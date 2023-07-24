# %% Imports
# Purpose: Functions for reading and writing COAPS data
from os.path import join
import xarray as xr
from datetime import datetime

# Only if debugging 
import sys
sys.path.append("../")  # (with interactive window)
# sys.path.append("eoas_pyutils")    # (without interactive window)


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


# %% Get HYCOM GoM data by date
def get_hycom_gom_raw_by_date(c_date, bbox=None):
    '''
    This code obtains the HYCOM GoM runs by date. Depending on the date selected is the run used to access the data. 
    This code only works if executed from the COAPS compute nodes (or if you have the folders mounted in your local machine)
    For 1993 - 2012 we use (50.1) https://www.hycom.org/data/gomu0pt04/expt-50pt1 /nexsan/archive/GOMu0.04_501/data/netcdf
    For 2013 we use (31.0) https://www.hycom.org/data/goml0pt04/expt-31pt0  /nexsan/archive/GOMl0.04_310
    For 2014 - 2018 (32.5) https://www.hycom.org/data/goml0pt04/expt-32pt5  /nexsan/archive/GOMl0.04_325
    For 2019 - Present (90.1m00) https://www.hycom.org/data/goml0pt04/expt-901m000  /nexsan/archive/GOMu0.04_901m000/data/hindcasts

    The names of the fiels are always the same:
    u for zonal velocity
    v for meridional velocity
    w for vertical velocity
    temperature for temperature
    depth for depth
    salinity for salinity
    ssh for surface elevation
    mld for mixed layer depth
    latitude for latitude
    longitude for longitude
    time for time
    '''

    # If the date is between 1993 and 2012
    year = c_date.year
    month = c_date.month
    day = c_date.day
    day_of_year = get_day_of_year_from_month_and_day(month, day, year)
    if  year >= 1993 and year <= 2012:
        print("Using HYCOM GoM 50.1 for years (1993 - 2012)")
        hycom_folder = '/nexsan/archive/GOMu0.04_501/data/netcdf'
        hycom_file_name = join(hycom_folder, str(year), f"hycom_gomu_501_{year}{month:02d}{c_date.day:02d}00_t*.nc")
        hycom_data = xr.open_mfdataset(hycom_file_name, decode_times=False)

        # Rename the variables
        hycom_data = hycom_data.rename({'water_u': 'u', 
                                        'water_v': 'v', 
                                        'water_temp': 'temperature', 
                                        'surf_el': 'ssh',
                                        'lat': 'latitude',
                                        'lon': 'longitude',
                                        })

    elif year == 2013:
        print("Using HYCOM GoM 31.0 for year 2013")
        hycom_folder = '/nexsan/archive/GOMl0.04_310'
        hycom_file_name = join(hycom_folder, str(year), f"hycom_gomu_501_{year}{month:02d}{c_date.day:02d}00_t000.nc")
    elif year >= 2014 and year <= 2018:
        print("Using HYCOM GoM 32.5 for years (2014 - 2018)")
        hycom_folder = '/nexsan/archive/GOMl0.04_325/data'
        day_of_year = get_day_of_year_from_month_and_day(month, c_date.day, year)
        hycom_file_name = join(hycom_folder, str(year), f"archv.{year}_{day_of_year:03d}_*_3z.nc")

        hycom_data = xr.open_mfdataset(hycom_file_name, decode_times=False)
        # Rename the variables
        hycom_data = hycom_data.rename({'MT': 'time',
                                        'Latitude': 'latitude',
                                        'Longitude': 'longitude',
                                        'w_velocity': 'w'
                                        })
    elif year >= 2019:
        print("Using HYCOM GoM 90.1m000 for years (2017 - 2023) hyndcast")
        hycom_folder = '/nexsan/archive/GOMu0.04_901m000/data/hindcasts'
        day_of_year = get_day_of_year_from_month_and_day(month, day, year)
        hycom_file_name = join(hycom_folder, str(year), f"hycom_gomu_901m000_{year}{month:02d}{day:02d}12_t*.nc")

        hycom_data = xr.open_mfdataset(hycom_file_name, decode_times=False)
        # Rename the variables
        hycom_data = hycom_data.rename({'water_u': 'u',
                                        'water_v': 'v',
                                        'water_temp': 'temperature',
                                        'surf_el': 'ssh',
                                        'lat': 'latitude',
                                        'lon': 'longitude',
                                        })


    # hycom_data = xr.load_dataset(hycom_file_name, decode_times=False)

    lats = hycom_data.latitude
    lons = hycom_data.longitude


    if bbox is not None:
        hycom_data = hycom_data.sel(latitude=slice(bbox[0],bbox[1]),
                                    longitude=slice(bbox[2], bbox[3])) 

    return hycom_data, lats, lons


# %% Main just for testing
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import numpy as np
    c_date = datetime(2019, 1, 10)
    c_date_str = c_date.strftime("%Y-%m-%d")
    bbox = [17.5, 32.5, -98, -76]

    hycom_data, lats, lons = get_hycom_gom_raw_by_date(c_date, bbox)

    # Make a grid 2x2 with the plots of surf_el, water_temp, salinity and velocity
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    # Superior title the date as a string
    fig.suptitle(c_date_str)
    # Plot surf_el
    axs[0, 0].imshow(hycom_data.ssh[0,:,:], origin="lower")
    axs[0, 0].set_title("Surface elevation")
    # Plot water_temp
    axs[0, 1].imshow(hycom_data.temperature[0,0,:,:], origin="lower")
    axs[0, 1].set_title("Water temperature")
    # Plot salinity in log scale
    axs[1, 0].imshow(hycom_data.salinity[0,0,:,:], origin="lower", norm=LogNorm())
    axs[1, 0].set_title("Salinity")
    # Plot velocity
    axs[1, 1].imshow(np.sqrt(hycom_data.u[0,0,:,:]**2 + hycom_data.v[0,0,:,:]**2), origin="lower")
    axs[1, 1].set_title("Velocity")
    plt.show()
    print("Done testing!")

# %%
