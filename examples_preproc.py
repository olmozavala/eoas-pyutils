import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import path
from shapely.geometry import Polygon

from viz_utils.eoa_viz import EOAImageVisualizer
from viz_utils.constants import PlotMode, BackgroundType
from proc_utils.comp_fields import vorticity, coriolis
from proc_utils.proj import haversineForGrid
from proc_utils.geometries import histogram_from_locations
from proc_utils.gom import lc_from_ssh

from shapely.geometry import LineString

## ============ Composite fields ===========
print("Reading data...")
input_file = "./test_data/hycom_gom.nc"
ds = xr.open_dataset(input_file, decode_times=False)
print(ds.info())
# Reading specific field and layers
lons = ds.lon
lats = ds.lat

## ------------ Vorticity -----------
vort = vorticity(ds.water_u, ds.water_v)
grid = np.meshgrid(lons, lats)
grid_dist = haversineForGrid(grid) / 1000
vort_norm = vorticity(ds.water_u, ds.water_v, grid_dist)
mag = np.sqrt(ds.water_u**2 + ds.water_v**2)
viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs", eoas_pyutils_path=".")
cbarlim = .3
viz_obj.plot_3d_data_npdict({'vort':vort[0,:]}, ['vort'], [0], 'Vorticity', 'myplot', mincbar=[-1*cbarlim], maxcbar=[cbarlim])
viz_obj.plot_3d_data_npdict({'mag':mag[0,:]}, ['mag'], [0], 'Magnitude', 'myplot')
viz_obj.plot_3d_data_npdict({'vortnorm':vort_norm[0,:]}, ['vortnorm'], [0], 'Vorticity normalized', 'myplot', mincbar=[-1*cbarlim], maxcbar=[cbarlim])

## ------------ Coriolis -----------
cor_par = coriolis(lats)
fig, ax = plt.subplots(1,1, figsize=(8,4))
ax.plot(lats, cor_par)
ax.set_title("Coriolis parameter")
ax.set_xlabel("Lats")
ax.set_ylabel("Coriolis")
plt.show()

## ============ Cropping fields ===========
ds_crop = ds.sel(lat=slice(24,30), lon=slice(-84, -78))  # Cropping by value
ds_crop = ds_crop.isel(time=0)  # Cropping by index
lats = ds_crop.lat
lons = ds_crop.lon
viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs", eoas_pyutils_path=".")
viz_obj.plot_3d_data_npdict(ds_crop, ['water_u'], [0], 'U', 'myplot', mincbar=[-.5], maxcbar=[.5])

## ============ Raster histogram ===========
file_name = "test_data/hycom_gom.nc"
ds = xr.open_dataset(file_name, decode_times=False)
lats = ds.lat.data
lons = ds.lon.data
viz_obj = EOAImageVisualizer(disp_images=True, output_folder='output', lats=lats, lons=lons)  # For Python console
viz_obj.plot_3d_data_npdict({'water_temp':ds.water_temp[0,:,:,:]}, ['water_temp'], title=F'Field Example (temp)', file_name_prefix='Test', z_levels=[0])

## ----------------- Histogram from locations
# Make grid of size ~2 degrees
gridres = 1
minlat, maxlat, rangelat = (18, 32, 32-18)
minlon, maxlon, rangelon = (-98, -76, 98-76)
lats_coarse = np.linspace(minlat, maxlat, int(rangelat/gridres)+1)
lons_coarse = np.linspace(minlon, maxlon, int(rangelon/gridres)+1)
ds_coarse = ds.interp(lat=lats_coarse, lon=lons_coarse) # Interpolate original temp grid to new resolution
grid_coarse = ds_coarse.water_temp[0, 0, :, :].data.copy() # Copy original grid
grid_coarse[~np.isnan(grid_coarse)] = 0
# Interpolate to new grid
test_lon = -93.0
test_lat = 22.0
# We make 6 locations and increase the number at center
# locations = zip([19, 20, 21, test_lat, test_lat, test_lat, ], [test_lon, test_lon, test_lon, -94, -95, -96])
locations = zip([19, 20, 21], [test_lon, test_lon, test_lon])
hist = histogram_from_locations(grid_coarse, lats_coarse, lons_coarse, locations)
viz_obj.plot_2d_data_np(hist, ['histogram'], title=F'Histogram locations', file_name_prefix='hist')

## ----------------- Intersection with geo_poly example ----------------------
geom_poly= np.array([[-87.5, 21.15], [-84.15, 22.35], [-82.9, 22.9], [-81, 22.9], [-81, 27], [-82.5, 32.5], [-76.5, 32.5], [-76.5, 16.5], [-90, 16.5], [-87.5, 21.15]])
polygon_shape = Polygon(geom_poly)

# viz_obj = EOAImageVisualizer(disp_images=True, output_folder='output', lats=lats, lons=lons, eoas_pyutils_path="../.", additional_polygons=[polygon_shape])
viz_obj.__setattr__('additional_polygons',[polygon_shape])
viz_obj.plot_3d_data_npdict({'water_temp':ds.water_temp[0,:,:,:]}, ['water_temp'], title=F'Field with geo_poly', file_name_prefix='Test', z_levels=[0])

print("Making the intersection...")
grid_bin = intersect_polygon_grid(ds.water_temp[0,0,:,:], lats, lons, geom_poly)
print("Done!")
viz_obj.plot_2d_data_np(grid_bin, ['binary_grid'], flip_data=False, rot_90=False, title=F'Intersection Example', file_name_prefix='Test')
print("Done")
##

