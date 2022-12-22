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
                grid[i, j] = gom_path.contains_points(([[lons[j], lats[i]]])).item()

    return grid

def histogram_from_locations(grid, lats, lons, locations):
    '''
    It computes the spatial histogram of a grid from a list of locations.
    Locations need to be specified (lat,lon) or (row, col). It assumes
    lons go from -180 to 180
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