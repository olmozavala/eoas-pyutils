# %% 
# INSTALL:
# mamba install conda-forge::copernicusmarine --yes
# copernicusmarine get --help

import os
# Get the current working directory
current_directory = os.getcwd()

# Check if 'download_data' directory exists in the current directory
if current_directory.find("download_data") != -1:
    os.chdir(current_directory.replace("download_data", ""))
else:
    print("'download_data' directory not found in the current working directory.")
print(f"Current working directory: {os.getcwd()}")

# %%
import copernicusmarine as cm
from download_data.Copernicus_Datasets import Copernicus_Datasets, Copernicus_Enum
import netrc
import datetime

print(f"Version of copernicusmarine: {cm.__version__}")

# %% Read credentials from .netrc file
secrets = netrc.netrc()
username, account, password = secrets.hosts['COPERNICUS']

def download_by_year(year, cop_ds, bbox, output_folder):
    start_date = datetime.date(year,1,1)
    end_date = datetime.date(year,12,31)
    cm.subset(
        dataset_id=cop_ds['id'],
        dataset_version=cop_ds['version'],
        variables= cop_ds['variables'],
        minimum_longitude=bbox[0],
        maximum_longitude=bbox[2],
        minimum_latitude=bbox[1],
        maximum_latitude=bbox[3],
        start_datetime=start_date.strftime("%Y-%m-%dT00:00:00"),
        end_datetime=end_date.strftime("%Y-%m-%dT00:00:00"),
        output_filename=f"{cop_ds['short_name']}_{start_date.year}.nc",
        output_directory=output_folder,
        username=username,
        password=password,
        force_download=True
    )


# %% -------- Download data by year ----------
bbox = (-99.0, 17.0, -74.0, 31.0) # This is the default BBOX for GoM
cop_ds = Copernicus_Datasets[Copernicus_Enum.CHLORA_L3_D]
for c_year in range(2022, 2024):
    # output_folder = "test_data/Satellite_Data_Examples/"
    # output_folder = "/unity/f1/ozavala/DATA/GOFFISH/CHLORA/COPERNICUS"
    output_folder = "~/Downloads"
    download_by_year(c_year, cop_ds, bbox, output_folder)

# %% TODO Understand: 
# export COPERNICUSMARINE_DISABLE_SSL_CONTEXT=True
# export COPERNICUSMARINE_MAX_CONCURRENT_REQUESTS=7
# %% 
# cm.subset(
#     dataset_id=cop_ds['id'],
#     dataset_version=cop_ds['version'],
#     variables= cop_ds['variables'],
#     minimum_longitude=bbox[0],
#     maximum_longitude=bbox[2],
#     minimum_latitude=bbox[1],
#     maximum_latitude=bbox[3],
#     start_datetime=start_date.strftime("%Y-%m-%dT00:00:00"),
#     end_datetime=end_date.strftime("%Y-%m-%dT00:00:00"),
#     output_filename=f"{cop_ds['short_name']}_{start_date.year}.nc",
#     output_directory=output_folder,
#     username=username,
#     password=password,
#     force_download=True
# )

