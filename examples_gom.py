import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import path
from shapely.geometry import Polygon
from shapely.geometry import LineString, Polygon

from viz_utils.eoa_viz import EOAImageVisualizer
from viz_utils.constants import PlotMode, BackgroundType
from proc_utils.comp_fields import vorticity, coriolis
from proc_utils.proj import haversineForGrid
from proc_utils.geometries import histogram_from_locations, intersect_polygon_grid
from proc_utils.gom import lc_from_ssh

## ============ Composite fields ===========
print("Reading data...")
input_file = "/home/olmozavala/Dropbox/MyProjects/OZ_LIB/eoas_preprocessing/Data/AVISO/SSH/2022-01.nc"
adt = xr.open_dataset(input_file, decode_times=False)
print(adt.info())
# Reading specific field and layers
lons = adt.longitude
lats = adt.latitude

##
from proc_utils.gom import lc_from_ssh
adt_slice = adt.adt[0,:,:]
lats_adt = adt.latitude
lons_adt = adt.longitude
lc = lc_from_ssh(adt_slice.values, lons_adt, lats_adt)

## Convert to list string of polygons

viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs",  show_var_names=True)
mylinestring = LineString(list(lc))
viz_obj.__setattr__('additional_polygons', [mylinestring])
viz_obj.plot_2d_data_np({'adt':adt}, ['adt'], 'LC', 'filepref', plot_mode=PlotMode.CONTOUR)
#
##

