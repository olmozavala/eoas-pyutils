from scipy.interpolate import interp1d
from skimage import measure
from matplotlib import path
import numpy as np

gom_bnd = np.array([[-87.5, 21.15], [-84.15, 22.35], [-82.9, 22.9], [-81, 22.9], [-81, 27], [-82.5, 32.5], [-76.5, 32.5], [-76.5, 16.5], [-90, 16.5], [-87.5, 21.15]])
gom_path = path.Path(gom_bnd)

def coriolis(lats):
    """ Computes the coriolis parameter from a list of lats """
    Omeg=2*np.pi/(24*3600);
    f=2*Omeg*np.sin(np.deg2rad(Y));
    return f


def change_units(cc, lon, lat):
    """skimage functions extracts the contour but the paths are image indices
    and not (lon, lat)"""
    flon = interp1d(np.arange(0,len(lon)), lon)
    flat = interp1d(np.arange(0,len(lat)), lat)
    newcc = []
    for cci in cc:
        try:
            t = np.zeros(cci.shape)
            t[:,0] = flat(cci[:,0])
            t[:,1] = flon(cci[:,1])
            newcc.append(t)
        except Exception as e:
            print(F"Failed for {cci} error: {e}")
    return newcc

#mean_adt=0.35869857046709513
def lc_from_ssh(ssh, lon, lat, mean_adt=0.35869857):
    """
    Obtains the Loop Current contour from ssh data
    https://duacs.cls.fr/faq/what-are-the-product-specification/different-sea-surface-heights-used-in-altimetry/
    Args:
        ssh:
        mean_adt:

    Returns:

    """
    ssh = ssh - mean_adt
    cco = measure.find_contours(ssh, 0.17, fully_connected='low', positive_orientation='low')
    cc = change_units(cco, lon, lat)
    # Remove contour that are eddies
    for i in range(len(cc) - 1, -1, -1):
        cci = cc[i]
        if np.linalg.norm(cci[0] - cci[-1]) < 0.5:  # presence of a loop
            del cc[i]

    #Remove contour that are outside of the Gom
    for i in range(len(cc) - 1, -1, -1):
        cci = cc[i]
        # print(cci.shape)
        if np.all(gom_path.contains_points(cci)):  # carribbean or atlantic ocean
            del cc[i]
        elif np.all(cci[:, 1] < -89):  # contour in the western gom
            del cc[i]
        elif len(cci) < 25:  # non-useful little contour
            del cc[i]

    indexes = [0]
    if len(cc) >= 2:
        indexes = [0, 1]
    lc_lats = np.concatenate([np.flip(cc[i][:,0]) for i in indexes])
    lc_lons = np.concatenate([np.flip(cc[i][:,1]) for i in indexes])
    pos = zip(lc_lons, lc_lats)
    return pos
