# %%
# Purpose: Functions for reading and writing COAPS data
from os.path import join
import numpy as np
import xarray as xr
from datetime import datetime, date, timedelta


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

def get_sss_by_date(sss_folder, c_date, bbox=None):
    '''
    Reads salinity single day for a given date. You can also specify a bounding box and the data will be cropped to that region.
    '''
    c_date_str = c_date.strftime("%Y%m%d")

    day_of_year = get_day_of_year_from_month_and_day(c_date.month, c_date.day, year=datetime.now().year)

    sss_file_name = join(sss_folder, str(c_date.year), f"RSS_smap_SSS_L3_8day_running_{c_date.year}_{day_of_year}_FNL_v05.0.nc")
    sss_data = xr.load_dataset(sss_file_name)
    if bbox is not None:
        sss_data = sss_data.sel( lat=slice(bbox[0],bbox[1]),
                                lon=slice(bbox[2],bbox[3]))

    lats = sss_data.lat
    lons = sss_data.lon

    return sss_data, lats, lons


# %%
# Testing the provided functions
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from io_utils.dates_utils import get_day_of_year_from_month_and_day

    aviso_folder = "/unity/f1/ozavala/DATA/GOFFISH/AVISO/GoM/"
    sst_folder = "/unity/f1/ozavala/DATA/GOFFISH/SST/OISST"
    c_date = datetime(2012, 1, 10)
    c_date_str = c_date.strftime("%Y-%m-%d")

    # Print 
    # aviso_data, lats, lons = get_aviso_by_month(aviso_folder, c_date, bbox=[17.5, 32.5, -98, -76])
    # print(f"Reading monthly aviso data... the number of times available are: {aviso_data.time.size} for date {c_date_str}")

    # aviso_data, lats, lons = get_aviso_by_date(aviso_folder, c_date, bbox=[17.5, 32.5, -98, -76])
    # print(f"Reading monthly aviso data... the number of times available are: {aviso_data.time.size} requested {c_date_str} available {aviso_data.time.values}")
    # plt.imshow(aviso_data.adt[:,:], origin="lower")
    # plt.show()

    # sst_data, lats, lons = get_sst_by_date(sst_folder, c_date, bbox=[17.5, 32.5, -98, -76])
    # print(f"Reading monthly aviso data... for date {c_date_str}")
    # plt.imshow(sst_data.analysed_sst[0,:,:], origin="lower")
    # plt.show()

    sss_data, lats, lons = get_sss_by_date(sst_folder, c_date, bbox=[17.5, 32.5, -98, -76])
    print(f"Reading monthly aviso data... for date {c_date_str}")
    plt.imshow(sss_data.analysed_sst[:,:], origin="lower")
    plt.show()