#!/usr/bin/env python
# coding: utf-8

# In[1]:


from __future__ import print_function
import os
import sys
import requests
import glob
import pyreadr
from pathlib import Path
import pandas as pd
import numpy as np
from io import StringIO
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib


# In[2]:


STARTYEAR = 2016
ENDYEAR   = 2023

NO_NAME = "NULL"
YEAR_LINE = "years:"
VAR_LINE  = "variables:"
DAYMET_VARIABLES = ['prcp']
DAYMET_YEARS     = [str(year) for year in range(STARTYEAR, ENDYEAR + 1)]
DAYMET_URL_STR = r'https://daymet.ornl.gov/single-pixel/api/data?lat={}&lon={}'
doy_files = glob.glob(str(Path("/home/phenology-urban/release/input/tree_list").expanduser() / f"*.rds"))


# In[3]:


def parse_params(line, param_list):
    start_idx = line.index(":") + 1
    line_split = line[start_idx:].split(",")
    requested_params = []
    for elem in line_split:
        if elem in param_list:
            requested_params.append(elem)
    return ",".join(requested_params)


# In[5]:


session = requests.Session()
MAX_WORKERS = 10

def daymet_url(lat, lon):
    return DAYMET_URL_STR.format(lat, lon) + var_str + years_str

def parse_and_aggregate(text_str, id_value):
    pos = text_str.find("year,yday,prcp")
    df = pd.read_csv(StringIO(text_str[pos:]), dtype={"year":"int32","yday":"int16"})
    prcp_col = next(c for c in df.columns if c.lower().startswith("prcp"))
    if prcp_col != "prcp": df = df.rename(columns={prcp_col: "prcp"})
    dt = pd.to_datetime(df["year"]*1000 + df["yday"], format="%Y%j")
    m  = dt.dt.month.values
    season = np.where(np.isin(m,[3,4,5]), "spring",
              np.where(np.isin(m,[6,7,8]), "summer",
              np.where(np.isin(m,[9,10,11]), "fall", "winter")))
    df["season"] = season
    df.loc[m == 12, "year"] += 1
    df["id"] = id_value
    out = (df.groupby(["id","year","season"], as_index=False)
             .agg(Sum_Sum_mm=("prcp","sum"),
                  prcp_count=("prcp", lambda s: s.notna().sum())))
    return out

def fetch_one(idx):
    lat, lon, name = lats[idx], lons[idx], names[idx]
    url = daymet_url(lat, lon)
    res = session.get(url, timeout=30)
    res.raise_for_status()
    return parse_and_aggregate(res.text, name)

requested_vars  = ",".join(DAYMET_VARIABLES)
requested_years = ",".join(DAYMET_YEARS)
var_str   = f"&measuredParams={requested_vars}" if requested_vars else ""
years_str = f"&year={requested_years}" if requested_years else ""

for doy_file in doy_files:
    data_list = pyreadr.read_r(doy_file)[None]
    city = os.path.splitext(os.path.basename(doy_file))[0]
    print(f"Processing city: {city}")

    out_dir  = Path.home() / "phenology-urban/release/input/climate/Daymet" / city
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "prcp_monthly_2016_2024.csv"
    if out_file.exists():
        print("File exists, skip:", out_file)
        continue

    names = data_list["id"].astype(str).tolist()
    lats  = data_list["lat"].astype(float).tolist()
    lons  = data_list["lon"].astype(float).tolist()

    all_results = []
    num_loc_requested = len(lats)

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = {ex.submit(fetch_one, i): i for i in range(num_loc_requested)}
        for n, fut in enumerate(as_completed(futures), 1):
            i = futures[fut]
            if i % 100 == 0:
                print("Processing:", names[i])
            try:
                out = fut.result()
                all_results.append(out)
            except Exception as e:
                print(f"[WARN] id={names[i]} failed: {e}")

    print("Finished downloading", len(all_results), "locations.")
    final_df = pd.concat(all_results, ignore_index=True)
    final_df.to_csv(out_file, index=False)  #
    print("Saved to:", out_file)

