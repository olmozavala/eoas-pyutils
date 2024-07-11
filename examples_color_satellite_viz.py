# %% [markdown]
# This notebook shows examples of different satellite images and provides a summary for each of them. 
%load_ext autoreload
%autoreload 2

# %%

from viz_utils.eoa_viz import EOAImageVisualizer
from viz_utils.constants import PlotMode, BackgroundType
import xarray as xr
import numpy as np
from matplotlib.colors import LogNorm


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

viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="output_imgs", eoas_pyutils_path=".",
                             contourf=False, show_var_names=True, BackgroundType=BackgroundType.WHITE)
viz_obj.plot_3d_data_npdict(aviso_df, ['adt','ugos','vgos','sla'], title="AVISO data", z_levels=[0])

# %% Sea Surface Temperature 
# NASA EarthData SMAP GHRSST L4 analysis 0.01 2002-June to Present 
 
# Read GHRSST data
# file = './test_data/Satellite_Data_Examples/GHRSST/2016/20160101120000-CMC-L4_GHRSST-SSTfnd-CMC0.1deg-GLOB-v02.0-fv03.0.nc'
# ghrsst_df = xr.open_dataset(file)
# ghrsst_df = ghrsst_df.sel(lat=slice(bbox[1], bbox[3]), lon=slice(bbox[0], bbox[2]))
# lats = ghrsst_df.lat.values
# lons = ghrsst_df.lon.values
# print(ghrsst_df)

# viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="output_imgs", eoas_pyutils_path=".",
#                              contourf=False, show_var_names=True, BackgroundType=BackgroundType.WHITE, coastline=True, 
#                              land=True)
# viz_obj.plot_3d_data_npdict(ghrsst_df, ['analysed_sst'], title="GHRSST", z_levels=[0],file_name_prefix="GHRSST")


# %% Chlorophyll-a
# Copernicus Marine Service Ocean Colour Global L3 0.0417° x 0.0417° (2002 to present)
file = './test_data/Satellite_Data_Examples/Chlora/Ocean_Color_2024-04-10_2024-04-12.nc'
oc_df = xr.open_dataset(file)
oc_df = oc_df.sel(latitude=slice(bbox[1], bbox[3]), longitude=slice(bbox[0], bbox[2]))
lats = oc_df.latitude.values
lons = oc_df.longitude.values

# %%
viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="output_imgs", eoas_pyutils_path=".",
                                contourf=False, show_var_names=True, background=BackgroundType.WHITE, coastline=True, 
                                land=True)
viz_obj.plot_3d_data_npdict(oc_df, ['CHL'], title="Copernicus Ocean Color ()", z_levels=[0], file_name_prefix="Ocean_Color",
                            norm=LogNorm(0.04, 1))

# %%
