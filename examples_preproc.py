import xarray as xr
import numpy as np

from viz_utils.eoa_viz import EOAImageVisualizer
from viz_utils.constants import PlotMode, BackgroundType
from proc_utils.comp_fields import vorticity
from proc_utils.proj import haversineForGrid
from proc_utils.gom import lc_from_ssh
from shapely.geometry import LineString

## ------------ Composite fields -----------
print("Reading data...")
input_file = "/home/olmozavala/Dropbox/TestData/netCDF/GoM/hycom_gomu_501_1993010100_t000.nc"
ds = xr.open_dataset(input_file, decode_times=False)
print(ds.info())
# Reading specific field and layers
lons = ds.lon
lats = ds.lat

# ## ------------ Vorticity -----------
# vort = vorticity(ds.water_u, ds.water_v)
# grid = np.meshgrid(lons, lats)
# grid_dist = haversineForGrid(grid) / 1000
# vort_norm = vorticity(ds.water_u, ds.water_v, grid_dist)
# mag = np.sqrt(ds.water_u**2 + ds.water_v**2)
# viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs", eoas_pyutils_path=".")
# cbarlim = .3
# viz_obj.plot_3d_data_npdict({'vort':vort[0,:]}, ['vort'], [0], 'Vorticity', 'myplot', mincbar=[-1*cbarlim], maxcbar=[cbarlim])
# viz_obj.plot_3d_data_npdict({'mag':mag[0,:]}, ['mag'], [0], 'Magnitude', 'myplot')
# viz_obj.plot_3d_data_npdict({'vortnorm':vort_norm[0,:]}, ['vortnorm'], [0], 'Vorticity normalized', 'myplot', mincbar=[-1*cbarlim], maxcbar=[cbarlim])
#
# ## ------------ Cropping fields -----------
# ds_crop = ds.sel(lat=slice(24,30), lon=slice(-84, -78))  # Cropping by value
# ds_crop = ds_crop.isel(time=0)  # Cropping by index
# lats = ds_crop.lat
# lons = ds_crop.lon
# viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs", eoas_pyutils_path=".")
# viz_obj.plot_3d_data_npdict(ds_crop, ['water_u'], [0], 'U', 'myplot', mincbar=[-.5], maxcbar=[.5])

## ---------- LC from GoM SSH ------------
file_name = "/Net/work/ozavala/GOFFISH/AVISO/1994-01.nc"
# lc = lc_from_ssh(ds.surf_el[0,:,:], lons, lats)
viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs",
                             eoas_pyutils_path=".", background=BackgroundType.BLUE_MARBLE_LR)
# viz_obj.__setattr__('additional_polygons', [LineString(lc)])
viz_obj.plot_3d_data_npdict(ds, ['surf_el'], [0], 'SSH', 'myplot', mincbar=[-1], maxcbar=[1] )

## ---------- LC from AVISO ------------
file_name = "/Net/work/ozavala/GOFFISH/AVISO/1994-01.nc"
ds_aviso_full = xr.open_dataset(file_name)
ds_aviso = ds_aviso_full.sel(latitude=slice(18,32), longitude=slice(-98, -78))  # Cropping by value
lons_aviso = ds_aviso.longitude
lats_aviso = ds_aviso.latitude
lc = lc_from_ssh(ds_aviso.adt[0,:,:], lons_aviso, lats_aviso)
viz_obj = EOAImageVisualizer(lats=lats_aviso, lons=lons_aviso, disp_images=True, output_folder="outputs",
                             eoas_pyutils_path=".", background=BackgroundType.BLUE_MARBLE_LR)
viz_obj.__setattr__('additional_polygons', [LineString(lc)])
viz_obj.plot_3d_data_npdict(ds_aviso, ['adt'], [0], 'SSH', 'myplot', mincbar=[-1], maxcbar=[1] )

##

