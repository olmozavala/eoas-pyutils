# %%
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from os.path import join
import os
import cmocean.cm as cm
import cartopy.crs as ccrs
from datetime import datetime, timedelta, date
import xarray as xr
from shapely.geometry.linestring import LineString
import sys
# sys.path.append("/home/olmozavala/Dropbox/MyProjects/OZ_LIB/eoas_preprocessing/")
sys.path.append("/unity/f1/ozavala/CODE/lce_ml_detection/eoas_pyutils") # Just when running this file directly for testing
from proc_utils.gom import lc_from_ssh

## Define dates and domain

# root_folder = "/Net/work/ozavala/GOFFISH/"
root_folder = "/unity/f1/ozavala/DATA/GOFFISH"

sss_path = "SSS/SMAP_Global"
adt_path = "AVISO"
sst_path = "SST/OISST"
chlora_path = "CHLORA/NOAA"

start_date = datetime.strptime("2020-01-01", "%Y-%m-%d")
end_date = datetime.strptime("2021-01-01", "%Y-%m-%d")
cur_date = start_date

# bbox = (-99, -74, 17, 31)
bbox= (-91, -80, 20, 28)
lon = (bbox[0], bbox[1])
lon_360 = (360 + bbox[0], 360 + bbox[1])
lat = (bbox[2], bbox[3])

# %%
def read_day(cur_date):
    '''
    Read and crop data for a given day
    Args:
        cur_date:
    Returns:
    '''
    year = cur_date.year
    month = cur_date.month
    day_year = cur_date.timetuple().tm_yday
    day_month = cur_date.day

    sss_file = join(root_folder, sss_path, str(year),
                    f"RSS_smap_SSS_L3_8day_running_{year}_{day_year:03d}_FNL_v05.0.nc")
    adt_file = join(root_folder, adt_path, f"{year}-{month:02d}.nc")
    sst_file = join(root_folder, sst_path, str(year),
                    f"{year}{month:02d}{day_month:02d}090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1_subset.nc")
    chlora_file = join(root_folder, chlora_path, f"{year}-{month:02d}-{day_month:02d}.nc")

    assert os.path.exists(sss_file)
    assert os.path.exists(adt_file)
    assert os.path.exists(sst_file)
    assert os.path.exists(chlora_file)

    # Reading data
    print("Reading...", flush=True)
    sss_orig = xr.load_dataset(sss_file)
    adt_orig = xr.load_dataset(adt_file)
    sst_orig = xr.load_dataset(sst_file)
    chlora_orig = xr.load_dataset(chlora_file)

    # Cropping to the same domain
    print("Cropping...", flush=True)
    adt = adt_orig.sel(latitude=slice(lat[0], lat[1]), longitude=slice(lon[0], lon[1]))
    sss = sss_orig.sel(lat=slice(lat[0], lat[1]), lon=slice(lon_360[0], lon_360[1]))
    sst = sst_orig.sel(lat=slice(lat[0], lat[1]), lon=slice(lon[0], lon[1]))
    chlora = chlora_orig.sel(latitude=slice(lat[1], lat[0]),
                             longitude=slice(lon[0], lon[1]))  # For some reason latitude is flipped

    sss_orig.close()
    adt_orig.close()
    sst_orig.close()
    chlora_orig.close()

    adt_slice = adt.adt[day_month-1, :, :]
    sss_slice = sss.sss_smap_40km[:, :]
    sst_slice = sst.analysed_sst[0, :, :]
    chlora_slice = chlora.chlor_a[0, 0, :, :]

    # Computing LC
    print("Computing LC...", flush=True)
    lats_adt = adt.latitude
    lons_adt = adt.longitude
    lc = lc_from_ssh(adt_slice.values, lons_adt, lats_adt)

    lons_lats = {'adt': (lons_adt, lats_adt), 'sss': (sss.lon, sss.lat), 'sst': (sst.lon, sst.lat), 'chlora': (chlora.longitude, chlora.latitude)}

    return adt_slice, sss_slice, sst_slice, chlora_slice, LineString(list(lc)), lons_lats

def plotData(fig, axs, adt, sss, sst, chlora, lc, lons_lats):
    print("Plotting...", flush=True)

    # Plot each slice imshow
    # im_adt = axs[0][0].imshow(adt, origin='lower', cmap=cm.balance, extent=bbox, transform=ccrs.PlateCarree())
    # im_sss = axs[0][1].imshow(sss, origin='lower', cmap=cm.haline, norm=LogNorm(vmin=30, vmax=38), extent=bbox, transform=ccrs.PlateCarree())
    # im_sst = axs[1][0].imshow(sst, origin='lower', cmap=cm.thermal, extent=bbox, transform=ccrs.PlateCarree())
    # im_chlor_a = axs[1][1].imshow(chlora, origin='upper', cmap=cm.algae, norm=LogNorm(vmin=.15, vmax=.3), extent=bbox, transform=ccrs.PlateCarree())

    # Plot each slice contourf
    lons, lats = lons_lats['adt']
    im_adt = axs[0][0].contourf(lons, lats, adt, 128, cmap=cm.balance)
    axs[0][0].contour(lons, lats, adt, color='b')
    lons, lats = lons_lats['sss']
    im_sss = axs[0][1].contourf(lons, lats, sss, 128, cmap=cm.haline, vmin=30, vmax=38)
    lons, lats = lons_lats['sst']
    im_sst = axs[1][0].contourf(lons, lats, sst, 128, cmap=cm.thermal)
    lons, lats = lons_lats['chlora']
    lev_exp = np.linspace(np.log10(.10), np.log10(.35), 128)
    levs = np.power(10, lev_exp)
    # print(levs)
    im_chlor_a = axs[1][1].contourf(lons, lats, chlora, 128, cmap=cm.algae, levels=levs)

    # Add LC
    x, y = lc.xy
    axs[0][0].plot(x, y, transform=ccrs.PlateCarree(), c='r', markersize=1)
    axs[0][1].plot(x, y, transform=ccrs.PlateCarree(), c='r', markersize=1)
    axs[1][0].plot(x, y, transform=ccrs.PlateCarree(), c='r', markersize=1)
    axs[1][1].plot(x, y, transform=ccrs.PlateCarree(), c='r', markersize=1)

    # Set titles
    axs[0][0].set_title("ADT")
    axs[0][1].set_title("SSS")
    axs[1][0].set_title("SST")
    axs[1][1].set_title("CHLORA ")

    # Assign gridlines
    for ax in fig.get_axes():
        gl = ax.gridlines(draw_labels=True, color='grey', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.left_labels = False

    # Add colorbars
    shrink = 0.7
    fig.colorbar(im_adt, ax=axs[0][0], shrink=shrink)
    fig.colorbar(im_sss, ax=axs[0][1], shrink=shrink)
    fig.colorbar(im_sst, ax=axs[1][0], shrink=shrink)
    fig.colorbar(im_chlor_a, ax=axs[1][1], shrink=shrink)

    plt.suptitle(cur_date)
    plt.tight_layout(pad=1.5)

    return fig

## Save plots as images
start_date = date(2019, 1, 1)
end_date = date(2020, 1, 1)

# cur_date = start_date
# while cur_date < end_date:
#     fig, axs = plt.subplots(2, 2, figsize=(20, 13), subplot_kw={'projection': ccrs.PlateCarree()})
#     file_name = join("/Net/work/ozavala/GOFFISH/imgs/SatelliteDataComparison", cur_date.strftime("%Y-%m-%d"))
#     adt, sss, sst, chlora, lc, lons_lats = read_day(cur_date)
#     fig = plotData(fig, axs, adt, sss, sst, chlora, lc, lons_lats)
#     plt.savefig(file_name)
#     # plt.show()
#     print(f"Done saving {cur_date}", flush=True)
#     cur_date += timedelta(days=1)
#     plt.close()
##

cur_date = date(2019, 5, 5)
adt, sss, sst, chlora, lc, lons_lats = read_day(cur_date)

##
fig, axs = plt.subplots(2, 2, figsize=(20, 13), subplot_kw={'projection': ccrs.PlateCarree()})
file_name = join("/Net/work/ozavala/GOFFISH/imgs/SatelliteDataComparison", cur_date.strftime("%Y-%m-%d"))
fig = plotData(fig, axs, adt, sss, sst, chlora, lc, lons_lats)
# plt.savefig(file_name)
plt.show()
print(f"Done saving {cur_date}", flush=True)
cur_date += timedelta(days=1)
plt.close()
##


# %% For years 1993 to 2023 plot the mean values for the SSH field
input_folder = join(root_folder, "AVISO/Global")

mean_ssh = []
for year in range(1993,2023):
    file = join(input_folder, f"{year}.nc")
    ds = xr.load_dataset(file)
    mean_ssh.append(ds.adt.mean())

# Plot
plt.plot(range(1993,2023), mean_ssh)
plt.show()
# %%
