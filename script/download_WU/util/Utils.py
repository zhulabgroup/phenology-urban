import config
from datetime import timedelta, date
import lxml.html as lh
import requests


class Utils:
    session = requests.Session()
    weather_station_url = None

    def __init__(self, session, weather_station_url):
        self.session = session
        self.weather_station_url = weather_station_url
        
    @classmethod
    def date_range_generator(cls, start, end = date.today()):
        current_date = start.replace(day=1)
        while current_date <= end:
            yield current_date
            if current_date.month == 12:
                current_date = current_date.replace(year=current_date.year + 1, month=1, day=1)
            else:
                current_date = current_date.replace(month=current_date.month + 1, day=1)

            

    @classmethod
    def date_url_generator(cls, weather_station_url, start_date, end_date = date.today()):
        date_range = Utils.date_range_generator(start_date, end_date)
        for date in date_range:
            date_string = date.strftime("%Y-%m-%d")
            url = f'{weather_station_url}/table/{date_string}/{date_string}/monthly'
            yield date_string, url
    
    @classmethod
    def date_url_array(cls, date_url_gen):
        date_url_arr = []
        for url in date_url_gen:
            date_url_arr.append(url)
        return date_url_arr

    @classmethod
    def fetch_data_table(cls, url):
        """
            Fetches a weather data url and checks if there are data entries for that date
        """
        html_string = cls.session.get(url)
        doc = lh.fromstring(html_string.content)
        
        data_table = doc.xpath('//*[@id="main-page-content"]/div/div/div/lib-history/div[2]/lib-history-table/div/div/div/table/tbody/tr')
        if data_table != []:
            return True
        else:
            return False

