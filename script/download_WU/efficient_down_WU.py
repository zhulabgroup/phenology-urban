#!/usr/bin/env python
# coding: utf-8

# In[1]:


import requests
import csv
import lxml.html as lh
import os
import time
import glob
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from os.path import expanduser

import config
from util.UnitConverter import ConvertToSystem
from util.Parser import Parser
from util.Utils import Utils
from pathlib import Path

# -------- 全局设置 --------
START_DATE = config.START_DATE
END_DATE = config.END_DATE
UNIT_SYSTEM = config.UNIT_SYSTEM
FIND_FIRST_DATE = config.FIND_FIRST_DATE

ADDRESS_DIR = Path(expanduser('~/phenology-urban/script/download_WU/address'))
station_files = ADDRESS_DIR.glob("*.txt")
MAX_WORKERS_PER_STATION = 10       
GLOBAL_MAX_CONNECTIONS = 64
city = Path(STATIONS_FILE).stem


# In[2]:


def make_session():
    s = requests.Session()
    retries = Retry(total=3, backoff_factor=0.3,
                    status_forcelist=(429, 500, 502, 503, 504))
    adapter = HTTPAdapter(max_retries=retries,
                          pool_connections=GLOBAL_MAX_CONNECTIONS,
                          pool_maxsize=GLOBAL_MAX_CONNECTIONS)
    s.mount('http://', adapter)
    s.mount('https://', adapter)
    s.headers.update({
        "User-Agent": "Mozilla/5.0 (compatible; WeatherScraper/1.0; +https://example.com/bot)"
    })
    return s

SESSION = make_session()

# 预编译 XPath（热点路径）
HISTORY_TBODY_XPATH = '//*[@id="main-page-content"]/div/div/div/lib-history/div[2]/lib-history-table/div/div/div/table/tbody'
CONDITION_CHECK_XPATH = '//*[@id="main-page-content"]/div/div/div/div[1]/div/span[2]/text()'
# lxml 用不了“预编译 XPath”API，这里只是把路径提取出来避免反复写字符串

FIELDNAMES = ['Date', 'Temp_High', 'Temp_Avg', 'Temp_Low',
              'DewPoint_High', 'DewPoint_Avg', 'DewPoint_Low',
              'Humidity_High', 'Humidity_Avg', 'Humidity_Low',
              'WindSpeed_High', 'WindSpeed_Avg', 'WindSpeed_Low',
              'Pressure_High', 'Pressure_Low', 'Precip_Sum']

def fetch_one_date(url, date_string, timeout=10):
    """抓取单个日期页面并解析为 rows（失败返回空列表，不抛异常阻断并发）。"""
    try:
        resp = SESSION.get(url, timeout=timeout)
        # 某些站点偶发空表，做 1-2 次轻量重试
        if resp.status_code >= 500:
            time.sleep(0.3)
            resp = SESSION.get(url, timeout=timeout)
        doc = lh.fromstring(resp.content)
        history_table = doc.xpath(HISTORY_TBODY_XPATH)
        if not history_table:
            return []  # 没有表格就跳过
        rows = Parser.parse_html_table(date_string, history_table)
        return rows or []
    except Exception:
        return []

def write_header_if_needed(csv_path):
    need_header = not os.path.exists(csv_path) or os.path.getsize(csv_path) == 0
    if need_header:
        with open(csv_path, 'a+', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
            if UNIT_SYSTEM == "metric":
                writer.writerow({
                    'Date': 'Date',
                    'Temp_High': 'HighTemp',
                    'Temp_Avg': 'AvgTemp',
                    'Temp_Low': 'LowTemp',
                    'DewPoint_High': 'HighDew',
                    'DewPoint_Avg': 'AvgDew',
                    'DewPoint_Low': 'LowDew',
                    'Humidity_High': 'High_%',
                    'Humidity_Avg': 'Avg_%',
                    'Humidity_Low': 'Low_%',
                    'WindSpeed_High': 'High_km/h',
                    'WindSpeed_Avg': 'Avg_km/h',
                    'WindSpeed_Low': 'Low_km/h',
                    'Pressure_High': 'High_hPa',
                    'Pressure_Low': 'Low_hPa',
                    'Precip_Sum': 'Sum_mm'
                })

def scrap_station(weather_station_url, city):
    station_name = weather_station_url.rstrip('/').split('/')[-1]
    out_dir  = Path.home() / "phenology-urban/data/raw/CONUS/WU" / city
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f'{station_name}.csv'

    # 已存在就跳过（保持你原逻辑）
    if os.path.exists(csv_path):
        print(f"{csv_path} exists, pass")
        return

    # 先检查站点状态页以避免 404/Offline 浪费
    try:
        r = SESSION.get(weather_station_url, timeout=5)
        if r.status_code == 404:
            print(f"404 Error detected for URL: {weather_station_url}")
            return
    except Exception as e:
        print(f"Pre-check error for {weather_station_url}: {e}")

    # 生成日期 URL
    url_gen = Utils.date_url_generator(weather_station_url, START_DATE, END_DATE)
    date_url_list = list(url_gen)  # [(date_string, url), ...]

    # 并发抓取该站点所有日期
    all_rows = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS_PER_STATION) as ex:
        futures = [ex.submit(fetch_one_date, url, date_string) for date_string, url in date_url_list]
        for fut in as_completed(futures):
            rows = fut.result()
            if rows:
                all_rows.extend(rows)

    # 写盘（一次性）
    if all_rows:
        write_header_if_needed(csv_path)
        with open(csv_path, 'a+', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
            writer.writerows(all_rows)
        print(f"[{station_name}] Saved {len(all_rows)} rows")


# In[ ]:


def main():
    ADDRESS_DIR = Path(expanduser('~/phenology-urban/script/download_WU/address'))
    txt_files = ADDRESS_DIR.glob("*.txt")
    for txt_file in txt_files[1:]:
        city = txt_file.stem
        print(f"=== Processing city: {city} ===")
        with open(txt_file, 'r', encoding='utf-8') as stations_file:
            for url in stations_file:
                url = url.strip()
                if not url:
                    continue
                print(url)
                scrap_station(url, city)

if __name__ == "__main__":
    main()

