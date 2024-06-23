# %% This file downloads data from NASA Ocean Color using wget
import os
from os.path import join
import datetime
from datetime import timedelta

start_date = datetime.datetime(2010, 1, 1)
end_date = datetime.datetime(2010, 1, 2)


output_folder = "./test_data/Satellite_Data_Examples/"

# EXAMPLE https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/JPSS1_VIIRS.20180121.L3m.DAY.CHL.chlor_a.4km.nc
# EXAMPLE https://oceandata.sci.gsfc.nasa.gov


for c_date in range(start_date, end_date):
    year = c_date.year
    output_folder = join(current_config.output_folder,str(year))
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    cmd = f"wget -r -l1 -H -t1 -nd -N -np -A.nc -erobots=off -i {c_date.strftime('%Y-%m-%d')}.txt"
    print(cmd)
    os.system(cmd)
    print(f"Output folder: {output_folder}")
    file_names = os.listdir(output_folder)
    for file_name in file_names:
        print(file_name)