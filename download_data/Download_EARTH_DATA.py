# %% This code is required at the beginning of the script to locate other modules
import os
# Get the current working directory
current_directory = os.getcwd()

# Check if 'download_data' directory exists in the current directory
if current_directory.find("download_data") != -1:
    print("here")
    os.chdir(current_directory.replace("download_data", ""))
else:
    print("'download_data' directory not found in the current working directory.")
print(f"Current working directory: {os.getcwd()}")

# %%
from os.path import join
from multiprocessing import Pool
import datetime
from datetime import timedelta
from io_utils.io_common import create_folder, dotdict, tuple_to_string
from subset_dataset_oisst import oisst
import subprocess

##% This program downloads data using PODAAC from EarthDAT
# 1) Install PODAAC pip install podaac-data-subscriber
# https://github.com/podaac/data-subscriber/blob/main/README.md



# %%
def par_download(proc_id):

    create_folder(current_config.output_folder)
    c_date = start_date
    c_end_date = c_date + timedelta(days=days_increment)
    i = 0
    while c_date < final_end_date:
        year = c_date.year
        output_folder = join(current_config.output_folder,str(year))
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        if i % TOT_PROC == proc_id:

            options = dotdict({'shortname': current_config["name"],
                       'start_date': c_date.strftime('%Y-%m-%dT%H:%M:%SZ'),
                       'end_date': c_end_date.strftime('%Y-%m-%dT%H:%M:%SZ'),
                       'box': bbox,
                       'onlySST': False,
                       'alltime': False})
# 
            cmd = f"podaac-data-downloader -c {options.shortname} -d {output_folder} --start-date {options.start_date} --end-date {options.end_date} -b=\"{tuple_to_string(bbox)}\""

            if subset:
                cmd += " --subset"

            print(cmd)
            os.system(cmd)
 
            file_names = os.listdir(output_folder)
            print(f"Output folder: {output_folder}")
            for file_name in file_names:
                try: 
                    if file_name.endswith(".nc") or file_name.endswith(".nc4"):
                        # TODO I'm not sure this renaming will work with all the files
                        new_name = current_config['rename_files'](file_name)
                        print("Renaming file from ", file_name, " to ", new_name)
                        os.rename(join(output_folder, file_name), join(output_folder, new_name))
                except Exception as e:
                    print(f"Error renaming file: {file_name}: {e}")

        c_date = c_date + timedelta(days=days_increment)
        c_end_date = c_date + timedelta(days=days_increment)
        i += 1

    print(F"Done all from process {proc_id}!")

# %%
if __name__ == '__main__':
    # Working configurations
    # def_output_folder = "/Net/work/ozavala/GOFFISH/"
    def_output_folder = "./test_data/Satellite_Data_Examples/"

    # ---------------------- GHRSST ----------------------
    # https://podaac.jpl.nasa.gov/dataset/MUR-JPL-L4-GLOB-v4.1  (analysis, 2002 to present)
    sst_ghrsst_v2_mur = dotdict({
        "name": "MUR-JPL-L4-GLOB-v4.1",
        "output_folder": join(def_output_folder,"SST", "GHRSST_MUR"),
        "rename_files": lambda file_name: sst_ghrsst_v2_mur["split_file"](file_name.split("_")[0].split("-")[0]),
        "split_file": lambda date_str: "GHRSST_MUR_" + date_str[:4] + "-" + date_str[4:6] + "-" + date_str[6:8] + ".nc"
    })
    # https://podaac.jpl.nasa.gov/dataset/CMC0.1deg-CMC-L4-GLOB-v3.0# (analysis, 2016 to present)
    sst_ghrsst_v4_cmc = dotdict({
        "name": "CMC0.1deg-CMC-L4-GLOB-v3.0",
        "output_folder": join(def_output_folder,"SST", "GHRSST_CMC"),
        "rename_files": lambda file_name: sst_ghrsst_v4_cmc["split_file"](file_name.split("_")[0].split("-")[0]),
        "split_file": lambda date_str: "GHRSST_CMC_" + date_str[:4] + "-" + date_str[4:6] + "-" + date_str[6:8] + ".nc"
    })

    # https://podaac.jpl.nasa.gov/dataset/SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5#capability-modal-download (2015 to present)
    sss_smap_v5 = dotdict({
        "name": "SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5",
        "output_folder": join(def_output_folder,"SSS", "SMAP_8day"),
        "rename_files": lambda file_name: sss_smap_v5["split_file"](file_name.split("_")[3]),
        "split_file": lambda date_str: "SSS_SMAP_v5_" + date_str[:4] + "-" + date_str[4:6] + "-" + date_str[6:8] + ".nc"
    })

    TOT_PROC = 1
    start_date = datetime.date(2016,1,1)
    final_end_date =  datetime.date(2016,1,2) # Testing single day
    # final_end_date = datetime.date.today()
    days_increment = 1
    bbox = (-99.0, 17.0, -74.0, 31.0) # This is the default BBOX for GoM
    subset = False  # If we want to crop to the BBOX while downloading
    # ---------- Execute single field
    current_config = sss_smap_v5 # It can be sst_ghrsst_v4_cmc or sss_smap_v5

    # Single day
    par_download(0)
    # Parallel
    # p = Pool(TOT_PROC)
    # p.map(par_download, range(TOT_PROC))

    # ---------- Execute all fields
    # for current_config in [sst_ghrsst_v2_mur, sst_ghrsst_v4_cmc, sss_smap_v5]:
    #     # Single day
    #     par_download(0)
    #     # Parallel
    #     # p = Pool(TOT_PROC)
    #     # p.map(par_download, range(TOT_PROC))
# %%
