import unicodedata
import re
from datetime import datetime


class Parser:
    @staticmethod
    def format_key(key: str) -> str:
        # Replace white space and delete dots
        return key.replace(' ', '_').replace('.', '')

    @staticmethod
    def parse_html_table(date_string: str, history_table: list) -> dict:

        table_rows = [tr for tr in history_table[0].xpath('//tr') if len(tr) == 16]
        headers_list = []
        data_rows = []

        column_mapping = {
            0: 'Date',
            1: 'Temp_High',
            2: 'Temp_Avg',
            3: 'Temp_Low',
            4: 'DewPoint_High',
            5: 'DewPoint_Avg',
            6: 'DewPoint_Low',
            7: 'Humidity_High',
            8: 'Humidity_Avg',
            9: 'Humidity_Low',
            10: 'WindSpeed_High',
            11: 'WindSpeed_Avg',
            12: 'WindSpeed_Low',
            13: 'Pressure_High',
            14: 'Pressure_Low',
            15: 'Precip_Sum'
        }

        # set Table Headers
        for header in table_rows[0]:
            headers_list.append(header.text)

        for tr in table_rows[1:]:

            row_dict = {}
            for i, td in enumerate(tr.getchildren()):
                td_content = re.sub(r"[°\xa0]", "",unicodedata.normalize("NFKD", td.text_content().strip()))
                column_name = column_mapping.get(i, f"Unknown_{i}")
                row_dict[column_name] = td_content

            data_rows.append(row_dict)
        
        return data_rows
