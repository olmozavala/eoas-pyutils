import os
from multiprocessing import Pool
import datetime
from datetime import timedelta
from pydap.client import open_url
import wget

## Downloads Chlorophyll Daily Merged VIIRS NPP and N20 9 km
# https://coastwatch.noaa.gov/cwn/cw_b/
# Not working

output_folder = "/Net/work/ozavala/GOFFISH/CHLORA/NOAA"
output_folder = "./DownloadedData/CHLORA/NOAA"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

lon = [262, 305]  # [-98, -60]
lat = [7.5, 50]

start_date = datetime.date(2010,2,10)
# final_end_date = datetime.date(2022,10,7)
final_end_date = datetime.date.today()
c_date = start_date

while c_date < final_end_date:
    date_str = c_date.strftime('%Y-%m-%d')
    print(f"Downloading file: {date_str}")
    url = f"https://coastwatch.noaa.gov//erddap/griddap/noaacwNPPN20S3ASCIDINEOFDaily.nc?chlor_a%5B({date_str}T12:00:00Z):1:({date_str}T12:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(50.0):1:(7.5)%5D%5B(-98):1:(-60)%5D"
    output_file = f"{output_folder}/{date_str}.nc"
    filename = wget.download(url, output_file)
    c_date += timedelta(days=1)