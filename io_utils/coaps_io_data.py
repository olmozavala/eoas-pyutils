# %% Imports
# Purpose: Functions for reading and writing COAPS data
from os.path import join
import numpy as np
import xarray as xr
from datetime import datetime

from io_utils.dates_utils import get_day_of_year_from_month_and_day

# %% AVISO by month
def get_aviso_by_month(aviso_folder, c_date, bbox=None):
    '''
    Reads AVISO monthly data for a given date. You can also specify a bounding box and the data will be cropped to that region.
    '''
    aviso_file_name = join(aviso_folder, f"{c_date.year}-{c_date.month:02d}.nc")
    aviso_data = xr.load_dataset(aviso_file_name)
    if bbox is not None:
        aviso_data = aviso_data.sel(latitude=slice(bbox[0],bbox[1]),
                                    longitude=slice(bbox[2],bbox[3]))

    lats = aviso_data.latitude
    lons = aviso_data.longitude

    return aviso_data, lats, lons

# %% AVISO by date
def get_aviso_by_date(aviso_folder, c_date, bbox=None):
    '''
    Reads AVISO single day for a given date. You can also specify a bounding box and the data will be cropped to that region.
    '''
    aviso_file_name = join(aviso_folder, f"{c_date.year}-{c_date.month:02d}.nc")
    aviso_data = xr.load_dataset(aviso_file_name)
    if bbox is not None:
        target_time = np.datetime64(f"{c_date.year}-{c_date.month:02d}-{c_date.day:02d}T00:00:00")
        # Calculate the time differences
        time_diff = np.abs(aviso_data["time"] - target_time)
        # Get the index of the closest time
        closest_index = np.argmin(time_diff.values)

        print("Closest time index:", closest_index)
        aviso_data = aviso_data.sel(time=aviso_data["time"][closest_index],
                                latitude=slice(bbox[0],bbox[1]),
                                longitude=slice(bbox[2],bbox[3]))

    lats = aviso_data.latitude
    lons = aviso_data.longitude

    return aviso_data, lats, lons

# %% SST by date
def get_sst_by_date(sst_folder, c_date, bbox=None):
    '''
    Reads SST single day for a given date. You can also specify a bounding box and the data will be cropped to that region.
    '''
    c_date_str = c_date.strftime("%Y%m%d")
    sst_file_name = join(sst_folder, str(c_date.year), f"{c_date_str}090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1_subset.nc")
    sst_data = xr.load_dataset(sst_file_name)
    if bbox is not None:
        sst_data = sst_data.sel( lat=slice(bbox[0],bbox[1]),
                                lon=slice(bbox[2],bbox[3]))

    lats = sst_data.lat
    lons = sst_data.lon

    return sst_data, lats, lons

# %% SSS by date
def get_sss_by_date(sss_folder, c_date, bbox=None):
    '''
    Reads salinity single day for a given date. You can also specify a bounding box and the data will be cropped to that region.
    '''
    c_date_str = c_date.strftime("%Y%m%d")

    day_of_year = get_day_of_year_from_month_and_day(c_date.month, c_date.day, year=datetime.now().year)

    sss_file_name = join(sss_folder, str(c_date.year), f"RSS_smap_SSS_L3_8day_running_{c_date.year}_{day_of_year:03d}_FNL_v05.0.nc")
    sss_data = xr.load_dataset(sss_file_name)
    if bbox is not None:
        sss_data = sss_data.sel( lat=slice(bbox[0],bbox[1]),
                                lon=slice((bbox[2] + 360)%360,(bbox[3] + 360)%360))

    lats = sss_data.lat
    lons = sss_data.lon

    return sss_data, lats, lons


# %% Chlora NOAA by date
def get_chlora_noaa_by_date(chlora_folder, c_date, bbox=None):
    '''
    Reads salinity single day for a given date. You can also specify a bounding box and the data will be cropped to that region.
    '''
    chlora_file_name = join(chlora_folder,  f"{c_date.year}-{c_date.month:02d}-{c_date.day:02d}.nc")
    chlora_data = xr.load_dataset(chlora_file_name)
    if bbox is not None:
        # TODO for some reason the latitude field is flipped
        chlora_data = chlora_data.sel( latitude=slice(bbox[1],bbox[0]),
                                       longitude=slice(bbox[2], bbox[3])) 


    lats = chlora_data.latitude
    lons = chlora_data.longitude

    return chlora_data, lats, lons