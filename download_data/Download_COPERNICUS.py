# %% 
# INSTALL:
# mamba install conda-forge::copernicusmarine
# copernicusmarine get --help

from calendar import monthrange
import os
import sys
sys.path.append("/unity/f1/ozavala/CODE/lce_ml_detection/eoas_pyutils") # Just when running this file directly for testing
# %%
import copernicusmarine as cm
# from download_data.Copernicus_Datasets import Copernicus_Datasets, Copernicus_Enum  # For running as a script
from Copernicus_Datasets import Copernicus_Datasets, Copernicus_Enum  # For Interactive window
import netrc
import datetime

print(f"Version of copernicusmarine: {cm.__version__}")

# %% Read credentials from .netrc file
secrets = netrc.netrc()
username, account, password = secrets.hosts['COPERNICUS']

# %% -------- Download data by year ----------
def download_by_year(year, cop_ds, bbox, output_folder):
    start_date = datetime.date(year,1,1)
    end_date = datetime.date(year,12,31)

    if cop_ds['short_name'] != '':
        output_filename = f"{cop_ds['short_name']}_{start_date.year}.nc"
    else:
        output_filename = f"{start_date.year}.nc"

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
    #     output_filename=output_filename,
    #     output_directory=output_folder,
    #     username=username,
    #     password=password,
    #     force_download=True
    # )

    # cm.subset(
    # dataset_id="cmems_obs-oc_glo_bgc-plankton_my_l3-olci-300m_P1D",
    # dataset_version="202211",
    # variables=["CHL", "CHL_uncertainty", "flags"],
    # minimum_longitude=bbox[0],
    # maximum_longitude=bbox[2],
    # minimum_latitude=bbox[1],
    # maximum_latitude=bbox[3],
    # start_datetime=start_date.strftime("%Y-%m-%dT00:00:00"),
    # end_datetime=end_date.strftime("%Y-%m-%dT00:00:00"),
    # force_download=True,
    # subset_method="strict",
    # disable_progress_bar=True,
    # )

    cm.subset(
  dataset_id="cmems_obs-oc_glo_bgc-plankton_my_l3-olci-300m_P1D",
  dataset_version="202211",
  variables=["CHL", "CHL_uncertainty", "flags"],
  minimum_longitude=-179.99722290039062,
  maximum_longitude=179.9972381591797,
  minimum_latitude=-89.99722290039062,
  maximum_latitude=89.99722290039062,
  start_datetime="2024-08-20T00:00:00",
  end_datetime="2024-08-20T00:00:00",
  force_download=True,
  subset_method="strict",
  disable_progress_bar=True,
)

# %% -------- Download data by month ----------
def download_by_month(year, cop_ds, bbox, output_folder):
    for month in range(1, 13):
        start_date = datetime.date(year,month,1)
        end_date = datetime.date(year, month, monthrange(year, month)[1])
        if cop_ds['short_name'] != '':
            output_filename = f"{cop_ds['short_name']}_{start_date.year}-{start_date.month:02d}.nc"
        else:
            output_filename = f"{start_date.year}-{start_date.month:02d}.nc"

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
            output_filename=output_filename,
            output_directory=output_folder,
            username=username,
            password=password,
            force_download=True
        )
        break

# %% -------- Download data by year ----------
bbox_gom = (-98.25, 7.25, -55.0, 50.0) # DO NOT DELETE THIS LINE, THIS BBOX ARE IMPORTANT!
bbox_global = (-180, -90, 180, 90) # DO NOT DELETE THIS LINE, THIS BBOX ARE IMPORTANT!

bbox = bbox_global

# cop_ds = Copernicus_Datasets[Copernicus_Enum.CHLORA_L3_D]
cop_ds = Copernicus_Datasets[Copernicus_Enum.SSH_DUACS_L4_D_1993]

# output_folder = "/unity/f1/ozavala/DATA/GOFFISH/CHLORA/COPERNICUS"
# output_folder = "/tmp/OZ/"
output_folder = "/unity/f1/ozavala/DATA/GOFFISH/CHLORA/COPERNICUS_L3"

for c_year in range(1997, 2024):
    download_by_year(c_year, cop_ds, bbox_global, output_folder)
    break
    # download_by_month(c_year, cop_ds, bbox, output_folder)

# %% TODO Understand: 
# export COPERNICUSMARINE_DISABLE_SSL_CONTEXT=True
# export COPERNICUSMARINE_MAX_CONCURRENT_REQUESTS=7