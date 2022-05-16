import numpy as np
from proc_utils.proj import haversine


def vorticity(u, v, dist_grid=None):
    '''
    It computes the vorticity between the u and v fields. If a distance grid is provided
    it considers it for the computation.
    :param u:
    :param v:
    :param dist_grid:
    :return:
    '''
    # ------- 2D ----------
    if len(u.shape) == 2:
        if dist_grid is None:
            vort = np.diff(v, axis=1)[:-1, :] - np.diff(u, axis=0)[:, :-1]
        else:
            vort = np.diff(v, axis=1)[:-1, :]/dist_grid[0] - np.diff(u, axis=0)[:, :-1]/dist_grid[1]

    # ------- 3D ----------
    if len(u.shape) == 3:  # Assumes first dimension is time
        if dist_grid is None:
            vort = np.diff(v, axis=2)[:,: -1, :] - np.diff(u, axis=1)[: , :, :-1]
        else:
            vort = np.diff(v, axis=1)[:,: -1, :] / dist_grid[0] - np.diff(u, axis=0)[:, :, :-1] / dist_grid[1]
    # ------- 4D ----------
    if len(u.shape) == 4:  # Assumes first two dimension are time and depth
        if dist_grid is None:
            vort = np.diff(v, axis=3)[:,:, : -1, :] - np.diff(u, axis=2)[:, :, :, :-1]
        else:
            vort = np.diff(v, axis=1)[:, :, : -1, :] / dist_grid[0] - np.diff(u, axis=0)[:, :, :, -1] / dist_grid[1]
    return vort
