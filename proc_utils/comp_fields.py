import numpy as np

def vorticity(u, v, grid_u=None, grid_v=None):
    # ------- 2D ----------
    if len(u.shape) == 2:
        if grid_u == None:
            vort = np.diff(v, axis=1)[:-1, :] - np.diff(u, axis=0)[:, :-1]
    # ------- 3D ----------
    if len(u.shape) == 3:  # Assumes first dimension is time
        if grid_u == None:
            vort = np.diff(v, axis=2)[:,: -1, :] - np.diff(u, axis=1)[: , :, :-1]
    # ------- 4D ----------
    if len(u.shape) == 4:  # Assumes first two dimension are time and depth
        if grid_u == None:
            vort = np.diff(v, axis=3)[:,:, : -1, :] - np.diff(u, axis=2)[:, :, :, :-1]
    return vort
