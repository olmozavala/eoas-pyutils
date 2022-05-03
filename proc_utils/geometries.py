from shapely.geometry import box, Polygon, Point, MultiPoint
from rtree import index

from matplotlib import path

def intersect_polygon_grid(lats, lons, polygon):
    '''
    It generates a binary grid with the intersected locations of a polygon
    Args:
        lats:
        lons:
        polygon:

    Returns:
    '''
    # probably can remove the loops here and do this one line contains_points and not masked!

    idx = index.Index()

    gom_path = path.Path(polygon)
    grid = np.full((len(lats),len(lons)), np.nan)
    # TODO it must be a better way to do it
    for i in range(0, len(lats)):
        for j in range(0, len(lons)):
            if not np.ma.is_masked(grid[i, j]):
                # grid[i, j] = polygon.intersection(Point(lons[i],lats[i])).area
                grid[i,j] = not gom_path.contains_points(([[lons[j], lats[i]]]))

    return grid

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
    # viz_obj.plot_3d_data_npdict({'water_temp':ds.water_temp[0,:,:,:]}, ['water_temp'], title=F'Field Example', file_name_prefix='Test', z_levels=[0])

    # ----------------- Intersection with polygon example ----------------------
    data = np.array([[-87.5, 21.15], [-84.15, 22.35], [-82.9, 22.9], [-81, 22.9], [-81, 27], [-82.5, 32.5], [-76.5, 32.5], [-76.5, 16.5], [-90, 16.5], [-87.5, 21.15]])
    polygon_shape = Polygon(data)

    viz_obj = EOAImageVisualizer(disp_images=True, output_folder='output', lats=lats, lons=lons, eoas_pyutils_path="../.", additional_polygons=[polygon_shape])
    # viz_obj.plot_3d_data_npdict({'water_temp':ds.water_temp[0,:,:,:]}, ['water_temp'], title=F'Field with polygon', file_name_prefix='Test', z_levels=[0])

    print("Making the intersection...")
    grid_bin = intersect_polygon_grid(lats, lons, data)
    print("Done!")
    viz_obj.plot_2d_data_np(grid_bin, ['binary_grid'], flip_data=False, rot_90=False, title=F'Intersection Example', file_name_prefix='Test')
    print("Done")