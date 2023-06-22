#%%
import xarray as xr
from shapely.geometry import LineString
from viz_utils.eoa_viz import EOAImageVisualizer
from viz_utils.constants import PlotMode


#%% ============ Composite fields ===========
print("Reading data...")
input_file = "/home/olmozavala/Dropbox/MyProjects/OZ_LIB/eoas_preprocessing/Data/AVISO/SSH/2022-02.nc"
adt = xr.open_dataset(input_file, decode_times=False)

bbox= (-91, -80, 20, 28)
lon = (bbox[0], bbox[1])
lon_360 = (360 + bbox[0], 360 + bbox[1])
lat = (bbox[2], bbox[3])
adt = adt.sel(latitude=slice(lat[0], lat[1]), longitude=slice(lon[0], lon[1]))

print(adt.info())
# Reading specific field and layers
lons = adt.longitude
lats = adt.latitude

#%%
from proc_utils.gom import lc_from_ssh
adt_slice = adt.adt[0,:,:]
lats_adt = adt.latitude
lons_adt = adt.longitude
lc = lc_from_ssh(adt_slice.values, lons_adt, lats_adt)

#%% Convert to list string of polygons

viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs",  show_var_names=True, eoas_pyutils_path=".")
mylinestring = LineString(list(lc))
viz_obj.__setattr__('additional_polygons', [mylinestring])
viz_obj.plot_2d_data_np(adt_slice, ['adt'], 'LC', 'filepref', plot_mode=PlotMode.RASTER)