#!/usr/bin/env python3

import argparse
import gzip
import threading
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from io import StringIO
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DAYMET_URL = "https://daymet.ornl.gov/single-pixel/api/data"
DAYMET_VARIABLES = ("prcp", "srad")
THREAD_LOCAL = threading.local()


def parse_daily_response(text, tree_id):
    lines = text.splitlines()
    header = next((index for index, line in enumerate(lines) if line.startswith("year,yday,")), None)
    if header is None:
        raise ValueError(f"Daymet response has no daily data header for tree {tree_id}")

    daily = pd.read_csv(StringIO("\n".join(lines[header:])), dtype={"year": "int32", "yday": "int16"})
    daily.columns = daily.columns.str.replace(r"\s*\(.*\)$", "", regex=True)
    required = {"year", "yday", *DAYMET_VARIABLES}
    if not required.issubset(daily.columns):
        raise ValueError(f"Daymet response is missing columns: {sorted(required - set(daily.columns))}")
    if daily.empty:
        raise ValueError(f"Daymet returned no daily data for tree {tree_id}")

    daily["date"] = pd.to_datetime(
        daily["year"].astype(str) + "-" + daily["yday"].astype(str).str.zfill(3),
        format="%Y-%j",
    ).dt.strftime("%Y-%m-%d")
    daily.insert(0, "id", tree_id)
    return daily.loc[:, ["id", "date", "year", "yday", *DAYMET_VARIABLES]]


def fetch_daily(tree_id, lat, lon, years):
    if not hasattr(THREAD_LOCAL, "session"):
        session = requests.Session()
        retries = Retry(total=3, backoff_factor=1, status_forcelist=(429, 500, 502, 503, 504))
        session.mount("https://", HTTPAdapter(max_retries=retries))
        THREAD_LOCAL.session = session
    response = THREAD_LOCAL.session.get(
        DAYMET_URL,
        params={
            "lat": lat,
            "lon": lon,
            "measuredParams": ",".join(DAYMET_VARIABLES),
            "year": ",".join(map(str, years)),
        },
        timeout=60,
    )
    response.raise_for_status()
    return parse_daily_response(response.text, tree_id)


def download_city(tree_file, output_dir, years, workers):
    import pyreadr

    city = tree_file.stem
    out_dir = output_dir / city
    out_file = out_dir / f"daymet_daily_{years[0]}_{years[-1]}.csv.gz"
    partial_file = out_file.with_name(out_file.name + ".partial")
    if out_file.exists():
        print(f"File exists, skipping {city}: {out_file}")
        return
    if partial_file.exists():
        raise FileExistsError(f"Incomplete file exists: {partial_file}; inspect it before retrying")

    trees = pyreadr.read_r(str(tree_file))[None]
    if not {"id", "lat", "lon"}.issubset(trees.columns):
        raise ValueError(f"{tree_file} must contain id, lat, and lon columns")
    trees = trees.loc[:, ["id", "lat", "lon"]].copy()
    trees["lat"] = pd.to_numeric(trees["lat"], errors="coerce")
    trees["lon"] = pd.to_numeric(trees["lon"], errors="coerce")
    valid = trees["id"].notna() & trees["lat"].between(-90, 90) & trees["lon"].between(-180, 180)
    if not valid.all():
        raise ValueError(f"{city} has {(~valid).sum()} trees with missing IDs or invalid coordinates")
    if trees["id"].duplicated().any():
        raise ValueError(f"{city} has duplicate tree IDs")
    if trees.empty:
        raise ValueError(f"{city} has no trees")
    trees["id"] = trees["id"].astype(str)

    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {len(trees)} trees for {city}")
    rows = iter(trees.itertuples(index=False, name=None))
    completed = 0
    with gzip.open(partial_file, "xt", encoding="utf-8", newline="") as output:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            pending = {}
            for tree_id, lat, lon in rows:
                pending[executor.submit(fetch_daily, tree_id, lat, lon, years)] = tree_id
                if len(pending) >= workers:
                    done, _ = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done:
                        tree_id = pending.pop(future)
                        try:
                            daily = future.result()
                        except Exception as error:
                            raise RuntimeError(f"Failed to download {city} tree {tree_id}") from error
                        daily.to_csv(output, index=False, header=completed == 0)
                        completed += 1
                        if completed % 100 == 0:
                            print(f"{city}: {completed}/{len(trees)} trees")
            for future, tree_id in pending.items():
                try:
                    daily = future.result()
                except Exception as error:
                    raise RuntimeError(f"Failed to download {city} tree {tree_id}") from error
                daily.to_csv(output, index=False, header=completed == 0)
                completed += 1

    partial_file.rename(out_file)
    print(f"Saved {completed} trees to {out_file}")


def main():
    parser = argparse.ArgumentParser(description="Download daily Daymet precipitation and shortwave radiation.")
    parser.add_argument("--tree-dir", type=Path, default= Path("/home/jializ/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/"),)
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "phenology-urban/data/raw/CONUS/Daymet/all_city_trees")
    parser.add_argument("--cities", nargs="+", help="City names matching RDS file stems; defaults to all cities")
    parser.add_argument("--start-year", type=int, default=2016)
    parser.add_argument("--end-year", type=int, default=2023)
    parser.add_argument("--workers", type=int, default=10)
    args = parser.parse_args()
    if args.start_year > args.end_year or args.workers < 1:
        parser.error("start-year must be <= end-year and workers must be at least 1")
    if not args.tree_dir.is_dir():
        parser.error(f"tree directory does not exist: {args.tree_dir}")
    files = sorted(args.tree_dir.glob("*.rds"))
    if args.cities:
        files = [args.tree_dir / f"{city}.rds" for city in args.cities]
        missing = [str(path) for path in files if not path.is_file()]
        if missing:
            parser.error(f"tree files not found: {', '.join(missing)}")
    if not files:
        parser.error(f"no RDS tree files found in {args.tree_dir}")
    years = range(args.start_year, args.end_year + 1)
    for tree_file in files:
        download_city(tree_file, args.output_dir, years, args.workers)


if __name__ == "__main__":
    main()
