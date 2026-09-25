# code/source_data/download_surveys.py
import os, sys
import json
import traceback
from datetime import datetime, timezone
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))
import argparse
from processing_help import retrieve_ehydro_data


def _dir_size_bytes(path: Path) -> int:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--districtSymbol", required=False)
    parser.add_argument("--channelArea", required=False)
    parser.add_argument("--surveyId", required=True)
    parser.add_argument("--startDate", required=True)
    parser.add_argument("--endDate", required=True)
    parser.add_argument("--outputDir", type=Path, required=True)
    parser.add_argument(
        "--maxWorkers",
        type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", 2)),
    )
    args = parser.parse_args()

    args.outputDir.mkdir(parents=True, exist_ok=True)

    done_file = args.outputDir / f"{args.surveyId}.done.json"
    failed_file = args.outputDir / f"{args.surveyId}.failed.json"

    # Idempotency: skip if already completed
    if done_file.exists():
        print(f"[skip] {args.surveyId} already downloaded (found {done_file})", flush=True)
        return

    # Clear any stale failure marker from a previous attempt
    if failed_file.exists():
        failed_file.unlink()

    try:
        surveys = retrieve_ehydro_data(
            data_dir=args.outputDir,
            district_symbol=args.districtSymbol,
            channel_area=args.channelArea,
            surveyID=args.surveyId,
            start_date=args.startDate,
            end_date=args.endDate,
            max_workers=args.maxWorkers,
        )
    except Exception as e:
        traceback.print_exc(file=sys.stderr)
        failed_file.write_text(json.dumps({
            "survey_id": args.surveyId,
            "status": "failed",
            "error": f"{type(e).__name__}: {e}",
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }, indent=2))
        sys.exit(1)

    done_file.write_text(json.dumps({
        "survey_id": args.surveyId,
        "status": "downloaded",
        "path": str(args.outputDir.resolve()),
        "size_bytes": _dir_size_bytes(args.outputDir),
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }, indent=2))
    print(f"[done] {args.surveyId} -> {args.outputDir}", flush=True)

    # NOTE: surveys would return a geodataframe with some unneeded properties of the returned surveys
    # leaving them in just because, but the retrieve_ehydro_data should be enough to download and unzip
    # the sounding data

if __name__ == "__main__":
    main()