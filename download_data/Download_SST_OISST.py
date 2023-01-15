import os
from os.path import join
from multiprocessing import Pool
import datetime
from datetime import timedelta
from io_utils.io_common import create_folder, dotdict
from subset_dataset_oisst import oisst
import subprocess


##% This program downloads SST data using PODAAC
# 1) Install PODAAC pip install podaac-data-subscriber
# https://github.com/podaac/data-subscriber/blob/main/README.md

# THE PODAAC CODE HAS BEEN UPDATED AND IT DOESN'T WORK ANYMORE. PLEASE READ THE DOCS AND UPDATE THIS FILES.
# ALSO MAKE A UNIT TEST TO  CATCH WHEN THINGS STOP WORKING

# You can't run this code from the Python Console (paths problems), run it 'normally'.

# root_output_folder = "/Net/work/ozavala/GOFFISH/SST/OISST"
root_output_folder = "./Data/SST/OISST"

TOT_PROC = 1
start_date = datetime.date(2011,12,1)
# final_end_date =  datetime.date(2011,12,2) # Testing single day
final_end_date = datetime.date.today()
days_increment = 1
bbox = (-99.0, -74.0, 17.0, 31.0) # This is the default BBOX for GoM

create_folder(root_output_folder)

def par_download(proc_id):
    c_date = start_date
    c_end_date = c_date + timedelta(days=days_increment)
    i = 0
    while c_date < final_end_date:
        year = c_date.year
        output_folder = join(root_output_folder,str(year))
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        if i % TOT_PROC == proc_id:
            options = dotdict({'shortname': 'MUR-JPL-L4-GLOB-v4.1',
                       'date0': c_date.strftime('%Y%m%d'),
                       'date1': c_end_date.strftime('%Y%m%d'),
                       'box': bbox,
                       'gridpoints': 1,
                       'onlySST': False,
                       'alltime': False})

            # cmd = F"./subset_dataset_oisst.py -s {c_date.strftime('%Y%m%d')} -f {c_end_date.strftime('%Y%m%d')} " \
            #       F"-b {bbox} -x MUR-JPL-L4-GLOB-v4.1"
            # print(cmd)
            # os.system(cmd)
            oisst(options)
            c_date = c_date + timedelta(days=days_increment)
            c_end_date = c_date + timedelta(days=days_increment)
            print("Done!")
        i += 1

        cmd = F"mv *.nc {output_folder}"
        print(cmd)
        os.system(cmd)

    print(F"Done all from process {proc_id}!")

## Single day
par_download(1)

## Parallel
p = Pool(TOT_PROC)
p.map(par_download, range(TOT_PROC))

