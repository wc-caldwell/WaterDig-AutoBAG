import os
from pygeohydro import EHydro
from concurrent.futures import ThreadPoolExecutor
from zipfile import ZipFile, BadZipFile
from pygeoogc.exceptions import ZeroMatchedError
import numpy as np
import verde as vd
import geopandas as gpd
import fiona
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors


def extract_single_zip(zipf):
    """Extract a single zip file and remove it after extraction"""
    try:
        extract_dir = zipf.parent / zipf.stem
        extract_dir.mkdir(parents=True, exist_ok=True)
        with ZipFile(zipf, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
        os.remove(zipf)
        return zipf.name, True
    except BadZipFile:
        os.remove(zipf)
        return zipf.name, False
    except Exception as e:
        return zipf.name, False

def retrieve_ehydro_data(data_dir, start_date, end_date, district_symbol=None, channel_area=None, surveyID = None, max_workers=8):
    """
    Function to retrieve eHydro survey data from the USACE database
    
    Parameters:
    -----------
    data_dir : Path
        Working directory for data storage
    start_date : str
        Start date in YYYY-MM-DD format
    end_date : str
        End date in YYYY-MM-DD format
    district_symbol : str
        USACE district code (e.g., 'CESWG')
    channel_area : str
        USACE NCF channel area named (e.g., 'CESWG_AN_01_BAY')
    surveyID : str
        USACE eHydro survey identifier (e.g., 'AN_01_BAY_20250523_CS')
    max_workers : int
        Number of parallel workers for zip extraction (default: 8)
    """

    where_parts = [
        f"surveydatestart >= '{start_date}'",
        f"surveydatestart <= '{end_date}'",
    ]

    id_only_parts = []

    if district_symbol:
        where_parts.append(f"usacedistrictcode = '{district_symbol}'")
        id_only_parts.append(f"usacedistrictcode = '{district_symbol}'")
        data_dir = data_dir / district_symbol

    survey_candidates = [None]
    if surveyID:
        parts = surveyID.split("_")
        # Right-to-left trimming: keep dropping trailing segments until a hit.
        survey_candidates = ["_".join(parts[:i]) for i in range(len(parts), 0, -1)]

    if channel_area:
        where_parts.append(f"channelareaidfk = '{channel_area}'")
        id_only_parts.append(f"channelareaidfk = '{channel_area}'")

    where_clause = " AND ".join(where_parts)
    id_only_clause = " AND ".join(id_only_parts) if id_only_parts else None

    def _query_with_optional_survey(base_parts, candidates, query_label):
        last_error = None
        attempted_queries = set()

        for candidate in candidates:
            survey_filters = [None]
            if candidate:
                # Try exact first, then prefix match for variable suffix tokens.
                safe_candidate = candidate.replace("'", "''")
                survey_filters = [
                    f"surveyjobidpk = '{safe_candidate}'",
                    f"surveyjobidpk LIKE '{safe_candidate}%'",
                ]

            for survey_filter in survey_filters:
                query_parts = list(base_parts)
                if survey_filter:
                    query_parts.append(survey_filter)
                query = " AND ".join(query_parts)

                if query in attempted_queries:
                    continue
                attempted_queries.add(query)

                print(f"{query_label}: {query}")
                try:
                    result = ehydro.bysql(query)
                    if candidate and candidate != surveyID:
                        print(f"Matched with trimmed SURVEY_ID: {candidate}")
                    if survey_filter and " LIKE " in survey_filter:
                        print("Matched via LIKE-prefix fallback on surveyjobidpk.")
                    return result
                except ZeroMatchedError as ex:
                    last_error = ex
                    continue

        if last_error is not None:
            raise last_error

        return ehydro.bysql(" AND ".join(base_parts))
    
    data_dir.mkdir(parents=True, exist_ok=True)
    
    ehydro = EHydro(data_type="outlines", cache_dir=data_dir)

    date_base_parts = [
        f"surveydatestart >= '{start_date}'",
        f"surveydatestart <= '{end_date}'",
    ]
    if district_symbol:
        date_base_parts.append(f"usacedistrictcode = '{district_symbol}'")
    if channel_area:
        date_base_parts.append(f"channelareaidfk = '{channel_area}'")

    id_base_parts = []
    if district_symbol:
        id_base_parts.append(f"usacedistrictcode = '{district_symbol}'")
    if channel_area:
        id_base_parts.append(f"channelareaidfk = '{channel_area}'")

    try:
        topobathy = _query_with_optional_survey(date_base_parts, survey_candidates, "Primary eHydro query")
    except ZeroMatchedError:
        # Date can be offset from the SURVEY_ID token date. Retry without date.
        if surveyID:
            print("No matches for date-constrained queries; retrying without date filters.")
            topobathy = _query_with_optional_survey(id_base_parts, survey_candidates, "Fallback eHydro query")
        else:
            raise

    topobathy.to_parquet(data_dir / 'ehydro.parquet')
    print(f'eHydro survey data saved locally to {data_dir}')
    print(f'Survey metadata saved to {data_dir / "ehydro.parquet"}')

    # Parallel zip extraction for speed
    zip_files = list(data_dir.glob('*.ZIP')) + list(data_dir.glob('*.zip'))
    if zip_files:
        print(f'Extracting {len(zip_files)} zip files in parallel...')
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(extract_single_zip, zip_files))
        extracted = sum(1 for _, success in results if success)
        print(f'Successfully extracted {extracted}/{len(zip_files)} files')

    return topobathy

def calculate_survey_characteristics(x, y, z):
    """
    Function to calculate the coefficient of variation of depth, density of the point cloud,
    and anisotropy in the distribution of the points within the point cloud. 
    """
    import numpy as np
    from shapely.geometry import MultiPoint
    
    # Density: total points / boundary area
    boundary = MultiPoint(list(zip(x, y))).convex_hull
    density = len(x) / boundary.area

    # Depth CV
    depth_cv = np.std(z, ddof=1) / abs(np.mean(z))

    # Anisotropy on full point cloud
    # this is anisotropy in the distribution of surveyed locations
    # for the anisotropy in the depth, directional variograms must be fit
    # this can be done quickly with GSTools using AIC, but that is a lot of work for 100 surveys
    cov = np.cov(x, y)
    eigs = np.sort(np.linalg.eigvalsh(cov))[::-1]
    anisotropy = (eigs[0] - eigs[1]) / (eigs[0] + eigs[1])

    return density, depth_cv, anisotropy

def decimate_survey_points(survey, tgt_num_pts, reduction_method='median', 
                           min_depth=0, max_depth=None, random_seed=1337):
    """
    Decimate survey points to exact target count.
    
    Parameters
    ----------
    survey_path : Path or str
        Path to the survey GDB file
    tgt_num_pts : int
        Target number of points after decimation
    reduction_method : str or callable, optional
        Method for block reduction ('median', 'mean', 'min', or 'nav_safe'), 
        default 'median'. Ignored if nav_safe=True.
    nav_safe : bool, optional
        If True, uses navigation-safe reduction (QC + shoal-biased minimum),
        default False
    min_depth : float, optional
        Minimum valid depth for QC filtering, default 0
    max_depth : float, optional
        Maximum valid depth for QC filtering, default None (no upper limit)
    random_seed : int, optional
        Random seed for reproducibility, default 1337
        
    Returns
    -------
    tuple
        (x, y, z, x_raw, y_raw, z_raw) arrays of decimated and raw coordinates
    """
    
    def _robust_minimum(values):
        """Shoal-biased reduction with outlier protection."""
        valid = values[np.isfinite(values)]
        if len(valid) == 0:
            return np.nan
        # 5th percentile protects against single erroneous shoal spikes
        # while still preserving real shoal features
        return np.percentile(valid, 5)
    
    def _find_optimal_spacing(x, y, z, target, reduction):
        """Binary search for optimal spacing."""
        low, high = 1.0, np.ptp(x) / 10
        
        while high - low > 0.1:
            mid = (low + high) / 2
            reducer = vd.BlockReduce(reduction=reduction, spacing=mid)
            (x_test, _), z_test = reducer.filter(coordinates=(x, y), data=z)
            
            if len(z_test) > target:
                low = mid
            else:
                high = mid
        
        return low
    
    # Load data
    if type(survey) == gpd.GeoDataFrame:
        # this is meant to let a user pass any geodataframe, as long as the xLocation, yLocation, and Z_use colums are defined in the geodataframe
        points = survey
    else:
        if 'SurveyPoint_HD' in fiona.listlayers(survey):
            points = gpd.read_file(survey, layer='SurveyPoint_HD')
        else:
            points = gpd.read_file(survey, layer='SurveyPoint')
        
    x_raw = points.xLocation.values
    y_raw = points.yLocation.values
    z_raw = points.Z_use.values
    
    # QC filtering (always applied if nav_safe, removes obvious errors)
    if reduction_method == 'nav_safe':
        valid_mask = np.isfinite(z_raw) & (z_raw >= min_depth)
        if max_depth is not None:
            valid_mask &= (z_raw <= max_depth)
        
        x_qc = x_raw[valid_mask]
        y_qc = y_raw[valid_mask]
        z_qc = z_raw[valid_mask]
        
        n_removed = len(z_raw) - len(z_qc)
        if n_removed > 0:
            print(f"QC removed {n_removed} invalid points ({100*n_removed/len(z_raw):.1f}%)")
        
        reduction = _robust_minimum
    else:
        x_qc, y_qc, z_qc = x_raw, y_raw, z_raw
        reduction = reduction_method
    
    # Find optimal spacing for block reduction
    spacing = _find_optimal_spacing(x_qc, y_qc, z_qc, tgt_num_pts, reduction)
    reducer = vd.BlockReduce(reduction=reduction, spacing=spacing)
    (x, y), z = reducer.filter(coordinates=(x_qc, y_qc), data=z_qc)
    
    # Remove duplicates
    coords_combined = np.column_stack([x, y])
    _, idx = np.unique(coords_combined, axis=0, return_index=True)
    x, y, z = x[idx], y[idx], z[idx]
    
    if len(z) > tgt_num_pts:
        # Calculate roughness for each decimated point
        tree = cKDTree(np.column_stack([x_qc, y_qc]))

        roughness = []
        for x_pt, y_pt in zip(x, y):
            indices = tree.query_ball_point([x_pt, y_pt], r=spacing*2)
            if len(indices) > 3:
                roughness.append(np.std(z_qc[indices]))
            else:
                roughness.append(0.1)

        roughness = np.array(roughness) + 0.01
        probs = roughness / roughness.sum()

        rng = np.random.default_rng(random_seed)
        idx = rng.choice(len(z), tgt_num_pts, replace=False, p=probs)
        x, y, z = x[idx], y[idx], z[idx]
    elif len(z) < tgt_num_pts:
        n_needed = tgt_num_pts - len(z)
        rng = np.random.default_rng(random_seed)
        # Sample from QC'd data, not raw, to maintain nav safety
        idx_extra = rng.choice(len(z_qc), n_needed, replace=False)
        x = np.concatenate([x, x_qc[idx_extra]])
        y = np.concatenate([y, y_qc[idx_extra]])
        z = np.concatenate([z, z_qc[idx_extra]])
    
    return x, y, z, x_raw, y_raw, z_raw

def validate_bathymetric_surface(surface_tif, x_raw, y_raw, z_raw, 
                                 output_csv=None, method_name=None):
    """
    Universal validation function - works with any bathymetric surface.
    
    Filters raw points to only those within raster extent before sampling.
    
    Parameters
    ----------
    surface_tif : str or Path
        Path to bathymetric surface GeoTIFF
    x_raw, y_raw, z_raw : array-like
        Raw survey point coordinates and depths
    output_csv : str or Path, optional
        Where to save results. Default: same location as TIF
    method_name : str, optional
        Method name for labeling (e.g., 'Median_TIN')
    
    Returns
    -------
    pandas.DataFrame
        Results with X, Y, True_Depth, Predicted_Depth, Residual
    dict
        Summary statistics
    """
    import pandas as pd
    import rasterio
    from pathlib import Path
    
    surface_tif = Path(surface_tif)
    
    print(f"\n{'='*60}")
    print(f"SURFACE VALIDATION")
    print(f"{'='*60}")
    print(f"Surface: {surface_tif.name}")
    if method_name:
        print(f"Method:  {method_name}")
    
    # Convert to numpy arrays
    x_raw = np.asarray(x_raw, dtype=np.float64)
    y_raw = np.asarray(y_raw, dtype=np.float64)
    z_raw = np.asarray(z_raw, dtype=np.float64)
    
    print(f"\nRaw Data:")
    print(f"  Total points:       {len(z_raw):,}")
    print(f"  X range:            {x_raw.min():.2f} to {x_raw.max():.2f}")
    print(f"  Y range:            {y_raw.min():.2f} to {y_raw.max():.2f}")
    print(f"  Z range:            {z_raw.min():.2f} to {z_raw.max():.2f} ft")
    
    # Open raster and get extent
    with rasterio.open(surface_tif) as src:
        raster_bounds = src.bounds
        raster_nodata = src.nodata
        raster_crs = src.crs
        
        print(f"\nRaster Properties:")
        print(f"  Bounds:             {raster_bounds}")
        print(f"  CRS:                {raster_crs}")
        print(f"  NoData value:       {raster_nodata}")
        print(f"  Shape:              {src.height} x {src.width}")
        print(f"  Resolution:         {src.res[0]:.2f} x {src.res[1]:.2f}")
        
        # Filter raw points to raster extent (with small buffer for edge cases)
        buffer = 0.1  # Small buffer to avoid edge exclusion
        extent_mask = (
            (x_raw >= raster_bounds.left - buffer) & 
            (x_raw <= raster_bounds.right + buffer) &
            (y_raw >= raster_bounds.bottom - buffer) & 
            (y_raw <= raster_bounds.top + buffer)
        )
        
        n_outside = (~extent_mask).sum()
        if n_outside > 0:
            print(f"\nFiltering:")
            print(f"  Points outside raster extent: {n_outside:,} ({100*n_outside/len(z_raw):.1f}%)")
        
        # Apply filter
        x_filtered = x_raw[extent_mask]
        y_filtered = y_raw[extent_mask]
        z_filtered = z_raw[extent_mask]
        
        print(f"  Points within extent:         {len(z_filtered):,}")
        
        if len(z_filtered) == 0:
            raise ValueError("No raw points overlap with raster extent!")
        
        # Sample surface at filtered raw locations
        coords = np.column_stack([x_filtered, y_filtered])
        preds = np.array([v[0] for v in src.sample(coords)], dtype=np.float32)
    
    # Check for nodata values
    if raster_nodata is not None:
        nodata_mask = (preds == raster_nodata)
        if nodata_mask.any():
            print(f"  Points sampled as NoData:     {nodata_mask.sum():,}")
            preds[nodata_mask] = np.nan
    
    # Calculate residuals
    resids = z_filtered - preds
    
    # Check for valid predictions
    valid_mask = ~np.isnan(resids) & ~np.isnan(preds)
    n_valid = valid_mask.sum()
    n_invalid = len(resids) - n_valid
    
    print(f"\nValidation Coverage:")
    print(f"  Valid predictions:            {n_valid:,} ({100*n_valid/len(resids):.1f}%)")
    print(f"  Invalid (NaN):                {n_invalid:,}")
    
    if n_valid == 0:
        raise ValueError("No valid predictions obtained!")
    
    if n_invalid > 0.1 * len(resids):
        print(f"  ⚠ Warning: {100*n_invalid/len(resids):.1f}% of points have invalid predictions")
    
    # Build results dataframe (full dataset including invalid)
    df = pd.DataFrame({
        'X': x_filtered,
        'Y': y_filtered,
        'True_Depth': z_filtered,
        'Predicted_Depth': preds,
        'Residual': resids
    })
    
    # Calculate statistics on valid data only
    valid_resids = resids[valid_mask]
    valid_preds = preds[valid_mask]
    valid_true = z_filtered[valid_mask]
    
    rmse = np.sqrt(np.mean(valid_resids**2))
    mae = np.mean(np.abs(valid_resids))
    bias = np.mean(valid_resids)
    std_dev = np.std(valid_resids)
    
    # Bootstrap CI for RMSE
    n_bootstrap = 1000
    subsample_size = min(10000, n_valid)
    rmse_samples = []
    for _ in range(n_bootstrap):
        sample = np.random.choice(valid_resids, size=subsample_size, replace=True)
        rmse_samples.append(np.sqrt(np.mean(sample**2)))
    rmse_ci = np.percentile(rmse_samples, [2.5, 97.5])
    rmse_se = np.std(rmse_samples)
    
    # Percentiles
    percentiles = {
        '5th': float(np.percentile(valid_resids, 5)),
        '25th': float(np.percentile(valid_resids, 25)),
        '50th': float(np.percentile(valid_resids, 50)),
        '75th': float(np.percentile(valid_resids, 75)),
        '95th': float(np.percentile(valid_resids, 95))
    }
    
    # Print statistics
    print(f"\n{'='*60}")
    print("VALIDATION STATISTICS")
    print(f"{'='*60}")
    print(f"\nError Metrics:")
    print(f"  RMSE:               {rmse:.3f} ± {rmse_se:.3f} ft")
    print(f"  RMSE 95% CI:        [{rmse_ci[0]:.3f}, {rmse_ci[1]:.3f}] ft")
    print(f"  MAE:                {mae:.3f} ft")
    print(f"  Bias:               {bias:.3f} ft")
    print(f"  Std Dev:            {std_dev:.3f} ft")
    print(f"\nResidual Distribution:")
    print(f"  Min:                {valid_resids.min():.3f} ft")
    print(f"  5th percentile:     {percentiles['5th']:.3f} ft")
    print(f"  25th percentile:    {percentiles['25th']:.3f} ft")
    print(f"  Median:             {percentiles['50th']:.3f} ft")
    print(f"  75th percentile:    {percentiles['75th']:.3f} ft")
    print(f"  95th percentile:    {percentiles['95th']:.3f} ft")
    print(f"  Max:                {valid_resids.max():.3f} ft")
    
    # Quick diagnostics
    print(f"\nQuick Diagnostics:")
    abs_resids = np.abs(valid_resids)
    print(f"  Points with |error| > 10 ft:  {(abs_resids > 10).sum():,} ({100*(abs_resids > 10).sum()/n_valid:.1f}%)")
    print(f"  Points with |error| > 20 ft:  {(abs_resids > 20).sum():,} ({100*(abs_resids > 20).sum()/n_valid:.1f}%)")
    print(f"  Points with |error| > 50 ft:  {(abs_resids > 50).sum():,} ({100*(abs_resids > 50).sum()/n_valid:.1f}%)")
    
    # Build metadata
    stats = {
        'method': method_name or surface_tif.stem,
        'surface_path': str(surface_tif),
        'n_raw_total': len(z_raw),
        'n_within_extent': len(z_filtered),
        'n_valid_predictions': int(n_valid),
        'n_invalid': int(n_invalid),
        'pct_valid': float(100 * n_valid / len(z_filtered)),
        'rmse': float(rmse),
        'rmse_se': float(rmse_se),
        'rmse_ci_95': [float(rmse_ci[0]), float(rmse_ci[1])],
        'mae': float(mae),
        'bias': float(bias),
        'std_dev': float(std_dev),
        'percentiles': percentiles,
        'residual_range': {
            'min': float(valid_resids.min()),
            'max': float(valid_resids.max())
        }
    }
    
    # Save results
    if output_csv is None:
        output_csv = surface_tif.parent / f"{surface_tif.stem}_validation.csv"
    else:
        output_csv = Path(output_csv)
    
    df.to_csv(output_csv, index=False)
    print(f"\nValidation results saved to: {output_csv}")
    print(f"{'='*60}\n")
    
    return df, stats