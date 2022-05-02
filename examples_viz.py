import sys
sys.path.append("hycom_utils/python")

from hycom.io import read_hycom_fields, subset_hycom_field, read_hycom_coords
from hycom.info import read_field_names
import xarray as xr
import numpy as np
from viz_utils.eoa_viz import EOAImageVisualizer

## ==================== 3D =========================================
## ------------ HYCOM -----------
# print("Reading data...")
# input_file = "./hycom_utils/test_data/archv.2009_153_00.a"
# coords_file = "./hycom_utils/test_data/regional.grid.a"
# layers = [0, 1, 2, 3]  # Depth layers we want to read
# print(F"The fields available are: {read_field_names(input_file)}")
# # Reading specific field and layers
# fields = ['srfhgt', 'temp', 'u-vel.']   # Which fields to plot
# hycom_fields = read_hycom_fields(input_file, fields, layers)
#
# print(F"The coords available are: {read_field_names(coords_file)}")
# coords = read_hycom_coords(coords_file, ['plon:', 'plat:'])
# lons = coords['plon']
# lats = coords['plat']
# print("Done!")
#
# viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs", eoas_pyutils_path=".")
# viz_obj.plot_3d_data_npdict(hycom_fields, ['temp','u-vel.'], [0], 'MyTitle', 'myplot')
# print("Done!")
#
# ## ------------ NetCDF -----------
# print("Reading data...")
# input_file = "/home/olmozavala/Dropbox/TestData/netCDF/ECMWF/ERA-interim/2012-08-01_2012-08-02.nc"
# df = xr.load_dataset(input_file)
# print(df.info())
# # Reading specific field and layers
# lons = df.longitude
# lats = df.latitude
#
# viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs", eoas_pyutils_path=".")
# viz_obj.plot_3d_data_npdict(df, ['u10','v10'], [0], 'MyTitle', 'myplot')
# print("Done!")

## ====================  2D =========================================
print("Reading data...")
input_file = "/home/olmozavala/Dropbox/TestData/netCDF/ECMWF/ERA-interim/2012-08-01_2012-08-02.nc"
df = xr.load_dataset(input_file)
print(df.info())
# Reading specific field and layers
lons = df.longitude
lats = df.latitude

viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs", eoas_pyutils_path=".")
npdata_2d = np.array([df.u10[0,:,:], df.v10[0,:,:]])
viz_obj.plot_2d_data_np(npdata_2d, ['u10','v10'], 'MyTitle', 'filepref')
print("Done!")
##

