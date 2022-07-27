from shapely.geometry import Polygon
import numpy as np

from matplotlib import path

def intersect_polygon_grid(grid, lats, lons, geo_poly):
    '''
    It generates a binary grid with the intersected locations of a geo_poly
    Args:
        lats:
        lons:
        geo_poly:

    Returns:
    '''
    # probably can remove the loops here and do this one line contains_points and not masked!
    gom_path = path.Path(geo_poly)
    # Pretty slow
    for i in range(0, len(lats)):
        for j in range(0, len(lons)):
            if not np.isnan(grid[i, j]):
                grid[i, j] = gom_path.contains_points(([[lons[j], lats[i]]]))

    return grid

def filter_locations_by_polygon(lats, lons, geo_poly):
    """
    Filters a list of coordinates by a polygon. Classic example, we need to know if
    the locations are inside our domain
    Args:
        lats:
        lons:
        geo_poly:

    Returns:
    """

if __name__ == "__main__":
    import xarray as xr
    import sys
    import numpy as np
    sys.path.append("eoas_pyutils/")
    from viz_utils.eoa_viz import EOAImageVisualizer

    file_name = "../test_data/hycom_gom.nc"
    ds = xr.open_dataset(file_name, decode_times=False)
    lats = ds.lat
    lons = ds.lon
    viz_obj = EOAImageVisualizer(disp_images=True, output_folder='output', lats=lats, lons=lons, eoas_pyutils_path="../.")
    viz_obj.plot_3d_data_npdict({'water_temp':ds.water_temp[0,:,:,:]}, ['water_temp'], title=F'Field Example', file_name_prefix='Test', z_levels=[0])

    ## ----------------- Intersection with geo_poly example ----------------------
    geom_poly= np.array([[-87.5, 21.15], [-84.15, 22.35], [-82.9, 22.9], [-81, 22.9], [-81, 27], [-82.5, 32.5], [-76.5, 32.5], [-76.5, 16.5], [-90, 16.5], [-87.5, 21.15]])
    polygon_shape = Polygon(geom_poly)

    viz_obj = EOAImageVisualizer(disp_images=True, output_folder='output', lats=lats, lons=lons, eoas_pyutils_path="../.", additional_polygons=[polygon_shape])
    viz_obj.plot_3d_data_npdict({'water_temp':ds.water_temp[0,:,:,:]}, ['water_temp'], title=F'Field with geo_poly', file_name_prefix='Test', z_levels=[0])

    print("Making the intersection...")
    grid_bin = intersect_polygon_grid(ds.water_temp[0,0,:,:], lats, lons, geom_poly)
    print("Done!")
    viz_obj.plot_2d_data_np(grid_bin, ['binary_grid'], flip_data=False, rot_90=False, title=F'Intersection Example', file_name_prefix='Test')
    print("Done")
##

