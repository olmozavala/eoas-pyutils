#%%
import xarray as xr
from shapely.geometry import LineString
from viz_utils.eoa_viz import EOAImageVisualizer
from viz_utils.constants import PlotMode
from proc_utils.gom import lc_from_ssh, lc_from_date

#%% ============ From specific AVISO file compute the LC
print("Reading data...")
input_file = "test_data/Satellite_data_Examples/aviso_2012-06.nc"
adt = xr.open_dataset(input_file, decode_times=False)

bbox= (-94, -80, 20, 30)
lon = (bbox[0], bbox[1])
lon_360 = (360 + bbox[0], 360 + bbox[1])
lat = (bbox[2], bbox[3])
adt = adt.sel(latitude=slice(lat[0], lat[1]), longitude=slice(lon[0], lon[1]))

print(adt.info())
# Reading specific field and layers
lons = adt.longitude
lats = adt.latitude

#%%
adt_slice = adt.adt[0,:,:]
lats_adt = adt.latitude
lons_adt = adt.longitude
lc = lc_from_ssh(adt_slice.values, lons_adt, lats_adt)

viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs",  show_var_names=True, eoas_pyutils_path=".")
mylinestring = LineString(list(lc))
viz_obj.__setattr__('additional_polygons', [mylinestring])
viz_obj.plot_2d_data_np(adt_slice, ['adt'], 'LC', 'filepref', plot_mode=PlotMode.RASTER)

#%% ============ From a date read the LC (from a pickle file inside COAPS)
import datetime
start_date = datetime.datetime(2001, 3, 1)

lc = lc_from_date(start_date)

viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs",  show_var_names=True, eoas_pyutils_path=".")
mylinestring = LineString(list(lc))
viz_obj.__setattr__('additional_polygons', [mylinestring])
viz_obj.plot_2d_data_np(adt_slice, ['adt'], 'LC', 'filepref', plot_mode=PlotMode.RASTER)

import datetime
start_date = datetime.datetime(2001, 3, 1)

lc = lc_from_date(start_date)

viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs",  show_var_names=True, eoas_pyutils_path=".")
mylinestring = LineString(list(lc))
viz_obj.__setattr__('additional_polygons', [mylinestring])
viz_obj.plot_2d_data_np(adt_slice, ['adt'], 'LC', 'filepref', plot_mode=PlotMode.RASTER)
# %%
