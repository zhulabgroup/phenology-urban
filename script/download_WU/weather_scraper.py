# Made with love by Karl
# Contact me on Telegram: @karlpy

import requests
import csv
import lxml.html as lh
import os
import config

from util.UnitConverter import ConvertToSystem
from util.Parser import Parser
from util.Utils import Utils

# configuration
stations_file = open('C:/Users/feixi/Desktop/the-weather-scraper-master copy/address/Allentown.txt', 'r')
URLS = stations_file.readlines()
# Date format: YYYY-MM-DD
START_DATE = config.START_DATE
END_DATE = config.END_DATE

# set to "metric" or "imperial"
UNIT_SYSTEM = config.UNIT_SYSTEM
# find the first data entry automatically
FIND_FIRST_DATE = config.FIND_FIRST_DATE


def scrap_station(weather_station_url):

    session = requests.Session()
    timeout = 10
    global START_DATE
    global END_DATE
    global UNIT_SYSTEM
    global FIND_FIRST_DATE
    
    url_gen = Utils.date_url_generator(weather_station_url, START_DATE, END_DATE)
    station_name = weather_station_url.split('/')[-1]
    file_name = f'{station_name}.csv'

    with open(file_name, 'a+', newline='') as csvfile:
        fieldnames = ['Date', 'Temp_High', 'Temp_Avg', 'Temp_Low',
                         'DewPoint_High', 'DewPoint_Avg', 'DewPoint_Low',
                         'Humidity_High', 'Humidity_Avg', 'Humidity_Low',
                         'WindSpeed_High', 'WindSpeed_Avg', 'WindSpeed_Low',
                         'Pressure_High', 'Pressure_Low', 'Precip_Sum']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # Write the correct headers to the CSV file
        if UNIT_SYSTEM == "metric":
            # 12:04 AM	24.4 C	18.3 C	69 %	SW	0.0 km/h	0.0 km/h	1,013.88 hPa	0.00 mm	0.00 mm	0	0 w/m²
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
        for date_string, url in url_gen:
            try:
                print(f'Scraping data from {url}')
                history_table = False
                while not history_table:
                    html_string = session.get(url, timeout=timeout)
                    doc = lh.fromstring(html_string.content)
                    history_table = doc.xpath('//*[@id="main-page-content"]/div/div/div/lib-history/div[2]/lib-history-table/div/div/div/table/tbody')
                    if not history_table:
                        print("refreshing session")
                        session = requests.Session()

                # parse html table rows
                data_rows = Parser.parse_html_table(date_string, history_table)

                if len(data_rows) != 0:
                    print(f'Saving {len(data_rows)} rows')
                writer.writerows(data_rows)
            except Exception as e:
                print(e)



for url in URLS:
    url = url.strip()
    station_name = url.split('/')[-1]
    file_name = f'{station_name}.csv'
    if os.path.exists(file_name):
        print(f"{file_name} exists, pass")
        continue
    session = requests.Session()
    html_station = session.get(url, timeout=5)
    if html_station.status_code == 404:
        print(f"404 Error detected for URL: {url}")
        continue
    condition = lh.fromstring(html_station.content)
    condition_check = condition.xpath('//*[@id="main-page-content"]/div/div/div/div[1]/div/span[2]/text()')
    print(condition_check[0])
#    if condition_check[0] == "Offline":
#        print(f"Offline site: {url}")
#        continue
    print(url)
    scrap_station(url)
