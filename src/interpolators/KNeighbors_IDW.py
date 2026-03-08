#!/usr/bin/env python
"""
IDW interpolation using scipy.spatial.cKDTree 
"""
import numpy as np
from pathlib import Path
from scipy.spatial import cKDTree
import geopandas as gpd
import matplotlib.pyplot as plt
import rasterio
from rasterio.transform import from_origin
from rasterio.features import geometry_mask
import gc
import multiprocessing
from joblib import Parallel, delayed


def _idw_chunk_worker(x_chunk, y_chunk, tree, z_data, n_neighbors, power):
    """
    Static function to interpolate a chunk using IDW.
    Defined at module level to avoid pickling issues with joblib.
    """
    xx, yy = np.meshgrid(x_chunk, y_chunk)
    grid_points = np.column_stack([xx.ravel(), yy.ravel()])
    
    # Query k nearest neighbors
    distances, indices = tree.query(grid_points, k=n_neighbors)
    
    # Handle coincident points (distance = 0)
    distances = np.maximum(distances, 1e-10)
    
    # Compute IDW weights: w = 1 / d^p
    weights = 1.0 / np.power(distances, power)
    weights_sum = weights.sum(axis=1, keepdims=True)
    weights_normalized = weights / weights_sum
    
    # Get neighbor values and compute weighted average
    neighbor_values = z_data[indices]
    interpolated = np.sum(weights_normalized * neighbor_values, axis=1)
    
    return interpolated.reshape(xx.shape).astype(np.float32)


class IDW:
    """Generate IDW bathymetric surfaces with optimized performance."""
    
    __slots__ = [
        'tgt_survey', 'survey_boundary', 'grid_res', 'tgt_num_pts',
        'reduction_method', 'plot_outputs', 'coordinates', 'bathymetry', 'x_grid',
        'y_grid', 'grid_z', 'transform', 'meta',
        'cv_metadata', 'output_path',
        '_crs', 'power', 'n_neighbors', 'n_jobs'
    ]
    
    def __init__(self, tgt_survey, reduction_method, x, y, z, grid_res=10.0, tgt_num_pts=5000, 
                 power=2, n_neighbors=32, plot_outputs=False, n_jobs=-1):
        
        self.tgt_survey = Path(tgt_survey)
        self.survey_boundary = gpd.read_file(tgt_survey, layer='SurveyJob')
        self._crs = self.survey_boundary.crs
        self.grid_res = grid_res
        self.tgt_num_pts = tgt_num_pts
        self.reduction_method = reduction_method
        self.power = power
        self.n_neighbors = n_neighbors
        self.plot_outputs = plot_outputs
        self.n_jobs = n_jobs if n_jobs != -1 else multiprocessing.cpu_count()
        
        # Store as float32 to reduce memory by 50%
        self.coordinates = (np.asarray(x, dtype=np.float32), np.asarray(y, dtype=np.float32))
        self.bathymetry = np.asarray(z, dtype=np.float32)

        self.x_grid = None
        self.y_grid = None
        self.grid_z = None
        self.transform = None
        self.meta = None
        self.cv_metadata = None
        self.output_path = None

    def _interpolate(self, chunk_size=500):
        """Perform IDW interpolation using scipy cKDTree."""
        x, y = self.coordinates
        
        # Pre-compute bounds
        x_min, x_max = float(x.min()), float(x.max())
        y_min, y_max = float(y.min()), float(y.max())
        x_range = x_max - x_min
        y_range = y_max - y_min
        
        # Build grid with consistent bounds (float32)
        nx = int(np.ceil(x_range / self.grid_res)) + 1
        ny = int(np.ceil(y_range / self.grid_res)) + 1
        self.x_grid = np.linspace(x_min, x_min + (nx - 1) * self.grid_res, nx, dtype=np.float32)
        self.y_grid = np.linspace(y_min, y_min + (ny - 1) * self.grid_res, ny, dtype=np.float32)
        
        # Build KD-tree from data points
        print(f"Building KD-tree with {len(self.bathymetry)} points...")
        data_points = np.column_stack([x, y])
        tree = cKDTree(data_points)
        
        # Initialize output grid
        self.grid_z = np.full((ny, nx), np.nan, dtype=np.float32)
        
        # Generate chunk jobs
        chunk_jobs = []
        for iy in range(0, ny, chunk_size):
            for ix in range(0, nx, chunk_size):
                iy_end = min(iy + chunk_size, ny)
                ix_end = min(ix + chunk_size, nx)
                x_chunk = self.x_grid[ix:ix_end]
                y_chunk = self.y_grid[::-1][iy:iy_end]  # [::-1] for correct orientation
                chunk_jobs.append((iy, ix, iy_end, ix_end, x_chunk, y_chunk))
        
        # Process chunks in parallel using loky backend
        # The worker function is at module level so it can be pickled
        print(f"Processing {len(chunk_jobs)} chunks using {self.n_jobs} cores...")
        results = Parallel(n_jobs=self.n_jobs, verbose=5, backend='loky')(
            delayed(_idw_chunk_worker)(
                x_chunk, y_chunk, tree, self.bathymetry, self.n_neighbors, self.power
            )
            for iy, ix, iy_end, ix_end, x_chunk, y_chunk in chunk_jobs
        )
        
        # Aggregate results
        for (iy, ix, iy_end, ix_end, _, _), z_chunk in zip(chunk_jobs, results):
            self.grid_z[iy:iy_end, ix:ix_end] = z_chunk
        
        # Transform
        self.transform = from_origin(
            x_min - self.grid_res / 2,
            y_max + self.grid_res / 2,
            self.grid_res, self.grid_res
        )
        
        # Apply boundary mask
        boundary_mask = geometry_mask(
            [self.survey_boundary.geometry.iloc[0]],
            out_shape=self.grid_z.shape,
            transform=self.transform,
            invert=True
        )
        self.grid_z[~boundary_mask] = np.nan
        del boundary_mask

    def _prepare_metadata(self):
        """Prepare COG metadata."""
        self.meta = {
            'driver': 'GTiff',
            'height': self.grid_z.shape[0],
            'width': self.grid_z.shape[1],
            'count': 1,
            'dtype': 'float32',
            'crs': self._crs,
            'transform': self.transform,
            'nodata': -9999.0,
            'tiled': True,
            'blockxsize': 512,
            'blockysize': 512,
            'compress': 'deflate',
            'predictor': 3,
            'zlevel': 6
        }

    def _write_output(self, output_path=None):
        """Write COG raster."""
        if output_path is None:
            output_path = self.tgt_survey.parent / f"IDW_{self.reduction_method}_{self.tgt_num_pts}.tif"
        else:
            output_path = Path(output_path)
        
        self.output_path = output_path
        
        with rasterio.open(output_path, 'w', **self.meta) as dst:
            dst.write(self.grid_z, 1)
            dst.build_overviews([2, 4, 8, 16, 32], rasterio.enums.Resampling.average)
        
        return output_path

    def generate(self, output_path=None):
        """Run complete IDW workflow."""
        self._interpolate()
        self._prepare_metadata()
        if self.plot_outputs == True:
            self.plot_output()
        out = self._write_output(output_path)
        
        # Aggressive cleanup after writing
        self.grid_z = None
        self.x_grid = None
        self.y_grid = None
        
        gc.collect()
        
        print(f"Exported: {out}")
        return out

    def plot_output(self):
        """Quick visualization."""
        if self.grid_z is None:
            raise ValueError("Grid data has been cleared. Call plot_output() before generate() completes.")
        fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
        extent = [float(self.x_grid.min()), float(self.x_grid.max()), 
                  float(self.y_grid.min()), float(self.y_grid.max())]
        im = ax.imshow(self.grid_z, origin="upper", extent=extent, cmap='viridis')
        plt.colorbar(im, ax=ax, label='Depth')
        ax.set(xlabel='Easting', ylabel='Northing', title='IDW Interpolation')
        ax.grid(color='white', linestyle='--', linewidth=0.5, alpha=0.3)
        plt.show()