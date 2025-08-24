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


# In[2]:


STARTYEAR = 2016
ENDYEAR   = 2023

NO_NAME = "NULL"
YEAR_LINE = "years:"
VAR_LINE  = "variables:"
DAYMET_VARIABLES = ['srad']
DAYMET_YEARS     = [str(year) for year in range(STARTYEAR, ENDYEAR + 1)]
DAYMET_URL_STR = r'https://daymet.ornl.gov/single-pixel/api/data?lat={}&lon={}'
doy_files = glob.glob(str(Path("/home/jializ/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input").expanduser() / f"*.rds"))


# In[3]:


def parse_params(line, param_list):
    start_idx = line.index(":") + 1
    line_split = line[start_idx:].split(",")
    requested_params = []
    for elem in line_split:
        if elem in param_list:
            requested_params.append(elem)
    return ",".join(requested_params)


# In[ ]:


for doy_file in doy_files:
    data_list = pyreadr.read_r(doy_file)[None]
    city = os.path.splitext(os.path.basename(doy_file))[0]
    print("Processing data for:", city)
    names = data_list["id"].astype(str).tolist()
    lats  = data_list["lat"].astype(float).tolist()
    lons  = data_list["lon"].astype(float).tolist()

    requested_vars = ",".join(DAYMET_VARIABLES)
    requested_years = ",".join(DAYMET_YEARS)

    var_str = ''
    if requested_vars:
        var_str = "&measuredParams=" + requested_vars

    years_str = ''
    if requested_years:
        years_str = "&year=" + requested_years
#         years_str = "&start=" + StartDate + "&end=" + EndDate
    
    out_dir = Path.home() / "phenology-urban/data/raw/CONUS/Daymet/all_city_trees" / city
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "srad_monthly_2016_2024.csv"
    if out_file.exists():
        print("File already exists, skipping:", out_file)
        continue

    num_loc_requested = len(lats)
    num_downloaded = 0
    all_results = []

    for i in range(0, num_loc_requested):
        curr_url = DAYMET_URL_STR.format(lats[i], lons[i]) + var_str + years_str
        if  i % 50 == 0:
            print("Processing:", names[i])
        res = requests.get(curr_url)
        if not res.ok:
            print("Could not access the following URL:", curr_url)
        else:
            text_str = res.content.decode("utf-8")
            lines = text_str.splitlines()
            for j, line in enumerate(lines):
                if line.startswith("year,yday,srad"):
                    start_idx = j
                    break

            clean_lines = lines[start_idx:]
            csv_str = "\n".join(clean_lines)
            df = pd.read_csv(StringIO(csv_str))
            df.columns = [c.strip() for c in df.columns]
            srad_col = next((c for c in df.columns if c.lower().startswith("srad")), None)
            if srad_col is None:
                raise ValueError("Not find srad column")
            df = df.rename(columns={srad_col: "srad"})
            df["year"] = df["year"].astype(int)
            df["yday"] = df["yday"].astype(int)
            df["date"] = pd.to_datetime(df["year"]*1000 + df["yday"], format="%Y%j")
            df["id"] = names[i]
            m = df["date"].dt.month
            df["season"] = np.select(
                [m.isin([3,4,5]), m.isin([6,7,8]), m.isin([9,10,11]), m.isin([12,1,2])],
                ["spring",        "summer",        "fall",           "winter"]
            )
            df.loc[m == 12, "year"] = df.loc[m == 12, "year"] + 1
            out = (df.groupby(["id", "year", "season"], as_index=False)
                     .agg(Avg_srad=("srad", "mean"),
                          srad_count=("srad", lambda s: s.notna().sum())))

            all_results.append(out)
            num_downloaded += 1

    print("Finished downloading", num_downloaded, "locations.")        
    final_df = pd.concat(all_results, ignore_index=True)
    final_df.to_csv(out_file, index=False)
    print("Saved to:", out_file)

