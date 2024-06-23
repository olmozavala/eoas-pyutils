# %% [markdown]
# This notebook shows examples of different satellite images and provides a summary for each of them. 

# %%
from viz_utils.eoa_viz import EOAImageVisualizer
from viz_utils.constants import PlotMode, BackgroundType
import xarray as xr
import numpy as np

# BBOX of the GoM
bbox = (-99.0, 17.0, -78.0, 31.0) # This is the default BBOX for GoM

# %% [markdown]
# # AVISO ADT
# AVISO ADT and Geostrophic velocities 0.25° × 0.25° (1993 to present)
# 
# The variables avaialble are:
# - Absolute Dynamic Topography (ADT)
# - Geostrophic velocities
# - Sea Level Anomaly (SLA)
# 
# https://resources.marine.copernicus.eu/product-detail/SEALEVEL_GLO_PHY_L4_MY_008_047/INFORMATION

# %%
file = './test_data/Satellite_Data_Examples/aviso_2012-06.nc'
aviso_df = xr.open_dataset(file)
# Crop to the Gulf of Mexico
aviso_df = aviso_df.sel(latitude=slice(bbox[1], bbox[3]), longitude=slice(bbox[0], bbox[2]))
# Print variables
print(aviso_df)

lats = aviso_df.latitude.values
lons = aviso_df.longitude.values

# viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="output_imgs", eoas_pyutils_path=".",
#                              contourf=False, show_var_names=True, BackgroundType=BackgroundType.WHITE)
# viz_obj.plot_3d_data_npdict(aviso_df, ['adt','ugos','vgos','sla'], title="AVISO data", z_levels=[0])

# %% [markdown]
# # Temperature 
# 
# NASA EarthData SMAP GHRSST L4 analysis 0.01 2002-June to Present 
# 
# https://podaac.jpl.nasa.gov/dataset/MUR-JPL-L4-GLOB-v4.1 
# 
# ## Install
# To use podac install with conda:
# 
# `conda install podaac-data-subscriber`
# 
# **You need** to have an account at https://urs.earthdata.nasa.gov/ and create a .netrc file in your home directory with your credentials. With format: 
# 
# ```
# machine urs.earthdata.nasa.gov
#     login <username>
#     password <password>
# ```

# %%
# Read GHRSST data
file = './test_data/Satellite_Data_Examples/GHRSST/2016/20160101120000-CMC-L4_GHRSST-SSTfnd-CMC0.1deg-GLOB-v02.0-fv03.0.nc'
ghrsst_df = xr.open_dataset(file)
ghrsst_df = ghrsst_df.sel(lat=slice(bbox[1], bbox[3]), lon=slice(bbox[0], bbox[2]))
lats = ghrsst_df.lat.values
lons = ghrsst_df.lon.values
print(ghrsst_df)

# Draw temperature data using directly cartopy, show the coastlines and the land
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature

# fig = plt.figure(figsize=(10, 5))
# ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
# ax.set_extent([bbox[0], bbox[2], bbox[1], bbox[3]], crs=ccrs.PlateCarree())
# ax.coastlines(resolution='10m')
# plt.pcolormesh(lons, lats, ghrsst_df.analysed_sst[0,:,:], transform=ccrs.PlateCarree())
# # Draw land on top
# ax.add_feature(cfeature.LAND, zorder=10, edgecolor='black')
# # Draw rivers and lakes
# plt.colorbar()
# plt.title("GHRSST")
# plt.show()

# %%
viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="output_imgs", eoas_pyutils_path=".",
                             contourf=False, show_var_names=True, BackgroundType=BackgroundType.WHITE, coastline=True, 
                             land=True)
viz_obj.plot_3d_data_npdict(ghrsst_df, ['analysed_sst'], title="GHRSST", z_levels=[0],file_name_prefix="GHRSST")
