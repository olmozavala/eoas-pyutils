import numpy as np
from proc_utils.proj import haversine

def coriolis(lats):
    """
    Obtains the coriolis parameter from a vector of latitudes
    Args:
        lats:

    Returns:

    """
    omeg = 2*np.pi/(24*3600)
    f = 2*omeg*np.sin(np.deg2rad(lats))
    return f


def vorticity(u, v, dist_grid=None):
    '''
    It computes the vorticity between the u and v fields. If a distance grid is provided
    it considers it for the computation.
    :param u:
    :param v:
    :param dist_grid:
    :return:
    '''
    final_vort = np.zeros(u.shape)
    # ------- 2D ----------
    if len(u.shape) == 2:
        if dist_grid is None:
            vort = np.diff(v, axis=1)[:-1, :] - np.diff(u, axis=0)[:, :-1]
        else:
            vort = np.diff(v, axis=1)[:-1, :]/dist_grid[0] - np.diff(u, axis=0)[:, :-1]/dist_grid[1]
        vort_dims = vort.shape
        final_vort[:vort_dims[0], :vort_dims[1]] = vort

    # ------- 3D ----------
    if len(u.shape) == 3:  # Assumes first dimension is time
        if dist_grid is None:
            vort = np.diff(v, axis=2)[:,: -1, :] - np.diff(u, axis=1)[: , :, :-1]
        else:
            vort = np.diff(v, axis=2)[:,: -1, :] / dist_grid[0] - np.diff(u, axis=1)[:, :, :-1] / dist_grid[1]

        vort_dims = vort.shape
        final_vort[:, :vort_dims[1], :vort_dims[2]] = vort

    # ------- 4D ----------
    if len(u.shape) == 4:  # Assumes first two dimension are time and depth
        if dist_grid is None:
            vort = np.diff(v, axis=3)[:,:, : -1, :] - np.diff(u, axis=2)[:, :, :, :-1]
        else:
            vort = np.diff(v, axis=3)[:, :, : -1, :] / dist_grid[0] - np.diff(u, axis=2)[:, :, :, :-1] / dist_grid[1]

        vort_dims = vort.shape
        final_vort[:, :, :vort_dims[2], :vort_dims[3]] = vort

    return final_vort