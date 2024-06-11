from enum import Enum

class PlotMode(Enum):
    """ Enum to select the different types of slices"""
    RASTER = 1   # Only plot the image (raster)
    CONTOUR = 2  # Only plot the contours
    MERGED = 3   # Plot raster and contours

class BackgroundType(Enum):
    """ Enum to select the different background options"""
    BLUE_MARBLE_LR = 1
    BLUE_MARBLE_HR = 2
    BATHYMETRY = 3
    TOPO = 4
    CARTO_DEF = 5
    BLACK = 6
    WHITE = 7
    GREY = 8
    NONE = 100

