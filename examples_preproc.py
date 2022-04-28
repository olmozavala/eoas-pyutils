import xarray as xr
import numpy as np

from viz_utils.eoa_viz import EOAImageVisualizer
from proc_utils.comp_fields import vorticity

## ------------ Composite fields -----------
print("Reading data...")
input_file = "/home/olmozavala/Dropbox/TestData/netCDF/GoM/hycom_gomu_501_1993010100_t000.nc"
ds = xr.open_dataset(input_file, decode_times=False)
print(ds.info())
# Reading specific field and layers
lons = ds.lon
lats = ds.lat

## ------------ Composite fields -----------
vort = vorticity(ds.water_u, ds.water_v)
mag = np.sqrt(ds.water_u**2 + ds.water_v**2)
viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs", eoas_pyutils_path=".")
viz_obj.plot_3d_data_npdict({'vort':vort[0,:]}, ['vort'], [0], 'Vorticity', 'myplot', mincbar=[-.5], maxcbar=[.5])
viz_obj.plot_3d_data_npdict({'mag':mag[0,:]}, ['mag'], [0], 'Magnitude', 'myplot')

## ------------ Cropping fields -----------
ds_crop = ds.sel(lat=slice(24,30), lon=slice(-84, -78))  # Cropping by value
ds_crop = ds_crop.isel(time=0)  # Cropping by index
lats = ds_crop.lat
lons = ds_crop.lon
viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs", eoas_pyutils_path=".")
viz_obj.plot_3d_data_npdict(ds_crop, ['water_u'], [0], 'U', 'myplot', mincbar=[-.5], maxcbar=[.5])

##

