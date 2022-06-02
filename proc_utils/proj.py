import numpy as np

def haversine(p1, p2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    All args must be of equal length.
    # Points in lat lon order
    """
    p1, p2 = map(np.radians, [p1, p2])

    dlat = p1[0] - p2[0]
    dlon = p1[1] - p2[1]

    a = np.sin(dlat/2.0)**2 + np.cos(p1[0]) * np.cos(p2[0]) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    dist = 6371000 * c # [m]
    return dist

def haversineForGrid(grid):
    """
    This function is used to obtain vertical and horizontal distances inside a grid
    :param grid:
    :return:
    """
    grid_rad = list(map(np.radians, grid))
    lat_rad = grid_rad[1]
    lon_rad = grid_rad[0]

    dlat = lat_rad[1:,:] - lat_rad[:-1,:]
    dlon = lon_rad[:,1:] - lon_rad[:,:-1]

    # We are creating horizontal and vertical distances
    out_dims = (2, grid[0].shape[0]-1,grid[0].shape[1]-1)
    output = np.zeros(out_dims)

    # Filling by cols
    for c_col in range(out_dims[2]):
        output[0,:,c_col] = np.sin(dlat[:,c_col]/2.0)**2
        output[1,:,c_col] = (np.cos(lat_rad[:,c_col+1]) * np.cos(lat_rad[:,c_col]) * np.sin(dlon[:,c_col]/2.0)**2)[:-1]

    c = 2 * np.arcsin(np.sqrt(output))
    dist = (6371000 * c)  # [m]
    return dist

if __name__ == "__main__":

    # ------ For two points
    # p1 = [0, 0]
    # p2 = [0, 1]
    # d = haversine(p1, p2)
    # print(F"Distance between {p1}(lat,lon) and {p2}(lat,lon) is {d:.3f} m")

    # ------ For grid
    lat = np.linspace(-10,10,21)
    lon = np.linspace(-20,20,41)
    grid = np.meshgrid(lon, lat)