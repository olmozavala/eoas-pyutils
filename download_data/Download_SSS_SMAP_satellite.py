# %%
import requests
import os
from os.path import join
from io_utils.io_common import create_folder

## This program dowloads SSS data. TODO the data is not cropped

def parallel_sss_download(output_folder, years, proc_id=0, TOT_PROC=1):
    for c_year in years:
        c_output_folder = join(output_folder,str(c_year))
        create_folder(c_output_folder)

        for c_day in range(1, 366):
            if c_day % TOT_PROC == proc_id:
                file_name = F"RSS_smap_SSS_L3_8day_running_{c_year}_{c_day:03d}_FNL_v05.0.nc"
                URL = F"https://data.remss.com/smap/SSS/V05.0/FINAL/L3/8day_running/{c_year}/{file_name}"
                try:
                    output_file = join(c_output_folder, file_name)
                    # ------- Only delete if the file is not the same size
                    if os.path.exists(output_file):
                        size = response.headers.get('content-length', 0)
                        size = int(size)
                        if os.path.getsize(output_file) == size:
                            continue
                        else:
                            os.remove(output_file)
                    print(F"Downloading file for day {c_year}-{c_day:03d}: {URL}")
                    response = requests.get(URL)

                    open(output_file, "wb").write(response.content)
                except Exception as e:
                    print(F"Failed for file: {file_name}")

    print("Done!")



# %%
