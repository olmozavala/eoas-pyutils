import sys
sys.path.append("hycom_utils/python")

from hycom.io import read_hycom_fields, subset_hycom_field, read_hycom_coords
from hycom.info import read_field_names
import xarray as xr
import numpy as np
from viz_utils.eoa_viz import EOAImageVisualizer
from viz_utils.constants import PlotMode
from shapely.geometry import LineString, Polygon
import cmocean.cm as ccm

## ==================== 3D =========================================
# ------------ HYCOM -----------
print("Reading data...")
input_file = "./hycom_utils/test_data/archv.2009_153_00.a"
coords_file = "./hycom_utils/test_data/regional.grid.a"
layers = [0, 1, 2, 3]  # Depth layers we want to read
print(F"The fields available are: {read_field_names(input_file)}")
# Reading specific field and layers
fields = ['srfhgt', 'temp', 'u-vel.']   # Which fields to plot
hycom_fields = read_hycom_fields(input_file, fields, layers)

##
print(F"The coords available are: {read_field_names(coords_file)}")
coords = read_hycom_coords(coords_file, ['plon:', 'plat:'])
lons = coords['plon']
lats = coords['plat']
print("Done!")

viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs", eoas_pyutils_path=".")

viz_obj.plot_3d_data_npdict(hycom_fields, ['temp','u-vel.'], [0], 'MyTitle', 'myplot')
hycom_fields['temp']
print("Done!")

## ------------ NetCDF 4D -----------
print("Reading data...")
input_file = "./test_data/hycom_gom.nc"
df = xr.load_dataset(input_file, decode_times=False)
print(df.info())
# Reading specific field and layers
lons = df.lon
lats = df.lat

viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs", eoas_pyutils_path=".")
viz_obj.plot_4d_data_npdict(df, ['water_temp', 'salinity'], z_levels=[0], title='Temp and Salinity', file_name_prefix='4d')
print("Done!")

## ====================  3D =========================================
print("Reading data...")
input_file = "./test_data/ecmwf_era_interim_2012-08-01.nc"
df = xr.load_dataset(input_file)
print(df.info())
# Reading specific field and layers
lons = df.longitude
lats = df.latitude

viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs", eoas_pyutils_path=".")
viz_obj.plot_3d_data_npdict(df, ['u10'], z_levels=[0], title='U', file_name_prefix='3d')
print("Done!")

## ====================  2D Changing plot mode =========================================
print("Reading data...")
input_file = "./test_data/ecmwf_era_interim_2012-08-01.nc"
df = xr.load_dataset(input_file)
print(df.info())
# Reading specific field and layers
lons = df.longitude
lats = df.latitude

viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs",
                             eoas_pyutils_path=".", show_var_names=True, contour_labels=[4])
npdata_2d = np.array([df.u10[0,:,:], df.v10[0,:,:]])

# ##
viz_obj.plot_2d_data_np(npdata_2d, ['u10','v10'], 'U and V', '2d_u_and_v')
viz_obj.plot_2d_data_np(npdata_2d, ['u10'], 'U Raster', '2d_raster', plot_mode=PlotMode.RASTER)
viz_obj.plot_2d_data_np(npdata_2d, ['u10'], 'U Contour', '2d_contour', plot_mode=PlotMode.CONTOUR)
viz_obj.plot_2d_data_np(npdata_2d, ['u10'], 'U Merged', '2d_merged', plot_mode=PlotMode.MERGED)
print("Done!")

## ========= Include vector shape ==============
# --- Linestring
mylinestring = LineString(((-80,10), (-82,20),(-90,20),(-93,30)))
viz_obj.__setattr__('additional_polygons', [mylinestring])
viz_obj.plot_2d_data_np(npdata_2d, ['u10'], 'With external linestring', 'filepref', plot_mode=PlotMode.CONTOUR)
## --- Polygon
mypolygon= Polygon(((-110,10),(-110,20),(-100,20),(-100,10),(-110,10)))
viz_obj.__setattr__('additional_polygons', [mypolygon])
viz_obj.plot_2d_data_np(npdata_2d, ['u10'], 'With external polygon', 'filepref', plot_mode=PlotMode.CONTOUR)

## ======== Stream plots ===========
x, y = np.meshgrid(lons, lats)
u = df.u10[0,:,:].data
v = df.v10[0,:,:].data
vel = np.sqrt(u**2 + v**2)
viz_obj.__setattr__('vector_field', {'u':u,'v':v,'x':x,'y':y})
viz_obj.plot_2d_data_np(npdata_2d, ['u10'], 'MyTitle', 'filepref', plot_mode=PlotMode.RASTER)

viz_obj.__setattr__('vector_field', {'u':u,'v':v,'x':x,'y':y, 'density':2, 'color':vel, 'cmap':ccm.phase})
viz_obj.plot_2d_data_np(npdata_2d, ['u10'], 'MyTitle', 'filepref', plot_mode=PlotMode.RASTER)

