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

def histogram_from_locations(grid, lats, lons, locations):
    '''
    It computes the spatial histogram of a grid from a list of locations.
    Locations need to be specified (lat,lon) or (row, col)
    :param grid:
    :param lats:
    :param lons:
    :param locations:
    :return:
    '''

    for row, col in locations:
        i = np.where(row <= lats)[0][0] - 1
        j = np.where(col <= lons)[0][0] - 1
        grid[i, j] += 1
    return grid

if __name__ == "__main__":
    import xarray as xr
    import sys
    import numpy as np
    sys.path.append("eoas_pyutils/")
    from viz_utils.eoa_viz import EOAImageVisualizer

    # file_name = "../test_data/hycom_gom.nc"
    file_name = "eoas_pyutils/test_data/hycom_gom.nc"
    ds = xr.open_dataset(file_name, decode_times=False)
    lats = ds.lat.data
    lons = ds.lon.data
    # viz_obj = EOAImageVisualizer(disp_images=True, output_folder='output', lats=lats, lons=lons, eoas_pyutils_path="../.")  # For F6
    viz_obj = EOAImageVisualizer(disp_images=True, output_folder='output', lats=lats, lons=lons, eoas_pyutils_path="eoas_pyutils")  # For Python console
    # viz_obj.plot_3d_data_npdict({'water_temp':ds.water_temp[0,:,:,:]}, ['water_temp'], title=F'Field Example', file_name_prefix='Test', z_levels=[0])

    # ----------------- Histogram from locations
    # Make grid of size ~2 degrees
    gridres = 1
    minlat, maxlat, rangelat = (18, 32, 32-18)
    minlon, maxlon, rangelon = (-98, -76, 98-76)
    lats_coarse = np.linspace(minlat, maxlat, int(rangelat/gridres)+1)
    lons_coarse = np.linspace(minlon, maxlon, int(rangelon/gridres)+1)
    ds_coarse = ds.interp(lat=lats_coarse, lon=lons_coarse)
    grid_coarse = ds_coarse.water_temp[0, 0, :, :].data.copy() # Copy original grid
    grid_coarse[~np.isnan(grid_coarse)] = 0
    # Interpolate to new grid
    # locations = zip([-93, -93.1, -95, -92, -93.1, -84], [25, 25, 26, 22, 25, 26])
    test_lon = -93.0
    test_lat = 22.0
    locations = zip([19, 20, 21, test_lat, test_lat, test_lat], [test_lon, test_lon, test_lon, -94, -95, -96])
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

