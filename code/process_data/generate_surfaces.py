# code/process_data/generate_surfaces.py
import os, sys
sys.path.insert(0, os.path.abspath('./src'))
import argparse
from pathlib import Path
import pandas as pd
import geopandas as gpd
import gstools as gs
import fiona
from datetime import datetime, timezone

from interpolators.TIN import TIN
from interpolators.NaturalNeighbor import NaturalNeighbor
from interpolators.KNeighbors_IDW import IDW
from interpolators.RBF import RBF
from interpolators.isoAutoKrige import isoAutoKrige

import json
import traceback
import argparse
import processing_help

#TODO: Add in argparse arguments. Will probably make options include the output directory, estimator, resolution, number of points to decimate to, and the cross validation approach

variograms = {
    "Gaussian": gs.Gaussian,
    "Exponential": gs.Exponential,
    "Matern": gs.Matern,
    "Stable": gs.Stable,
    "Rational": gs.Rational,
    "Circular": gs.Circular,
    "Spherical": gs.Spherical,
    "SuperSpherical": gs.SuperSpherical,
    "JBessel": gs.JBessel,
}

reduction_method = "median"

# TODO: Remove this from the interpolator classes
# since executing in HPC, cannot check the matplotlib outputs
plots = False


def _dir_size_bytes(path: Path) -> int:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variogramEstimator", type=str, required=True)
    parser.add_argument("--gridResolution", type=float, required=True)
    parser.add_argument("--soundings", type=Path, required=True)
    parser.add_argument("--chunkedKriging", required=True)
    parser.add_argument("--numPoints", type=int, required=True)
    parser.add_argument(
        "--maxWorkers",
        type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", 2)),
    )
    args = parser.parse_args()

    done_file = Path(args.soundings).parent / f"{Path(args.soundings).stem}.done.json"
    failed_file = Path(args.soundings).parent / f"{Path(args.soundings).stem}.failes.json"

    # Idempotency: skip if already completed
    if done_file.exists():
        print(f"[skip] {args.surveyId} already downloaded (found {done_file})", flush=True)
        return

    # Clear any stale failure marker from a previous attempt
    if failed_file.exists():
        failed_file.unlink()

    try:
        if 'SurveyPoint_HD' in fiona.listlayers(args.soundings):
            survey_pts = gpd.read_file(args.soundings, layer = 'SurveyPoint_HD')
        else:
            survey_pts = gpd.read_file(args.soundings, layer = 'SurveyPoint')

        x = survey_pts.xLocation.values
        y = survey_pts.yLocation.values
        z = survey_pts.Z_use.values

        TINinterpolator = TIN(
            tgt_survey = args.soundings,
            reduction_method = reduction_method,
            x=x, y=y, z=z,
            grid_res=args.gridResolution,
            tgt_num_pts=args.numPoints,
            plot_outputs=plots,
            n_jobs=args.maxWorkers
            )

        tin_path = TINinterpolator.generate()
        # _,_ = processing_help.validate_bathymetric_surface(tin_path, x_raw, y_raw, z_raw)

        NatNinterpolator = NaturalNeighbor(
                tgt_survey = args.soundings, 
                reduction_method = reduction_method,
                x=x, y=y, z=z, 
                grid_res=args.gridResolution,
                tgt_num_pts=args.numPoints,
                plot_outputs=plots,
                n_jobs=args.maxWorkers
            )

        natn_path = NatNinterpolator.generate()
        # _,_ = processing_help.validate_bathymetric_surface(natn_path, x_raw, y_raw, z_raw)

        IDWinterpolator = IDW(
                tgt_survey=args.soundings,
                reduction_method=reduction_method,
                x=x, y=y, z=z,
                grid_res=args.gridResolution,
                tgt_num_pts=args.numPoints,
                plot_outputs=plots,
                n_jobs=args.maxWorkers
            )

        idw_path = IDWinterpolator.generate()
        # _,_ = processing_help.validate_bathymetric_surface(idw_path, x_raw, y_raw, z_raw)

        RBFinterpolator = RBF(
                tgt_survey=args.soundings,
                reduction_method=reduction_method,
                x=x, y=y, z=z, 
                grid_res=args.gridResolution,
                tgt_num_pts=args.numPoints,
                plot_outputs=plots,
                n_jobs=args.maxWorkers
            )

        rbf_path = RBFinterpolator.generate(auto_tune=True)
        # _,_ = processing_help.validate_bathymetric_surface(rbf_path, x_raw, y_raw, z_raw)

        for detrend_bool in [False, True]:
        
            isoKrige_AIC = isoAutoKrige(
                                tgt_survey = args.soundings, 
                                vario_models=variograms,
                                reduction_method=reduction_method,
                                x=x, y=y, z=z, 
                                spline_damping=1e-1,                                            
                                vario_estimator = args.variogramEstimator,            
                                chunked_kriging = False,
                                grid_res=args.gridResolution,
                                detrend=detrend_bool,
                                tgt_num_pts=args.numPoints,
                                plot_outputs=plots,
                                n_jobs=args.maxWorkers
                        )
            krige_path_AIC, _ = isoKrige_AIC.generate()
            # _,_ = processing_help.validate_bathymetric_surface(krige_path_AIC, x_raw, y_raw, z_raw)

    except Exception as e:
        traceback.print_exc(file=sys.stderr)
        failed_file.write_text(json.dumps({
            "survey_id": Path(args.soundings).stem,
            "status": "failed",
            "error": f"{type(e).__name__}: {e}",
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }, indent=2))
        sys.exit(1)

    done_file.write_text(json.dumps({
        "survey_id": Path(args.soundings).stem,
        "status": "surfaces generated",
        "path": str(Path(args.soundings).parent),
        "size_bytes": _dir_size_bytes(Path(args.soundings).parent),
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }, indent=2))
    print(f"[done] {Path(args.soundings).stem} surfaces -> {Path(args.soundings).parent}", flush=True)

    # NOTE: surveys would return a geodataframe with some unneeded properties of the returned surveys
    # leaving them in just because, but the retrieve_ehydro_data should be enough to download and unzip
    # the sounding data

if __name__ == "__main__":
    main()