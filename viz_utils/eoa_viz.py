from os import listdir
from os.path import join

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm
from proj_utils.common import create_folder
from viz_utils.constants import PlotMode, BackgroundType
import pylab
import numpy as np
import cmocean

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy

def select_colormap(field_name):
    cmaps_fields = []
    if field_name in ('ssh', 'srfhgt'):
        # cmaps_fields.append(cmocean.cm.deep_r)
        cmaps_fields.append(cmocean.cm.curl)
    elif field_name in ('temp', 'sst', 'temperature'):
        cmaps_fields.append(cmocean.cm.thermal)
    elif field_name == "salin" or field_name == "sss" or field_name == "sal":
        cmaps_fields.append(cmocean.cm.haline)
    elif field_name in ('u', 'v', 'u-vel.', 'v-vel.'):
        cmaps_fields.append(cmocean.cm.speed)
    elif field_name.find("error") != -1:
        cmaps_fields.append(cmocean.cm.diff)
    return cmaps_fields


class EOAImageVisualizer:
    """This class makes plenty of plots assuming we are plotting Geospatial data (maps).
    It is made to read xarrays, numpy arrays, and numpy arrays in dictionaries"""
    _COLORS = ['y', 'r', 'c', 'b', 'g', 'w', 'k', 'y', 'r', 'c', 'b', 'g', 'w', 'k']
    _figsize = 8
    _font_size = 30
    _fig_prop = 1.8  # Proportion of each figure w/h
    _units = ''
    _max_imgs_per_row = 4
    _mincbar = np.nan  # User can set a min and max colorbar values to 'force' same color bar to all plots
    _maxcbar  = np.nan
    _flip_data = True
    _eoas_pyutils_path =  './eoas_pyutils'# This is the path where the eoas_utils folder is stored with respect to the main project
    _contourf = False  # When plotting non-regular grids and need precision
    _background = BackgroundType.BLUE_MARBLE_LR  # Select the background to use
    _auto_colormap = True  # Selects the colormap based on the name of the field
    _show_var_names = False  # Includes the name of the field name in the titles

    def __init__(self, disp_images=True, output_folder='output',
                 lats=[-90,90], lons =[-180,180],
                 projection=ccrs.PlateCarree(), **kwargs):
        # All the arguments that are passed to the constructor of the class MUST have its name on it.
        self._disp_images = disp_images
        self._output_folder = output_folder
        self._projection = projection
        self._extent = self.getExtent(lats, lons)
        self._lats = lats
        self._lons = lons
        for arg_name, arg_value in kwargs.items():
            self.__dict__["_" + arg_name] = arg_value

    def __getattr__(self, attr):
        '''Generic getter for all the properties of the class'''
        return self.__dict__["_" + attr]

    def __setattr__(self, attr, value):
        '''Generic setter for all the properties of the class'''
        self.__dict__["_" + attr] = value

    def add_colorbar(self, fig, im, ax, show_color_bar, label=""):
        # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.colorbar.html
        if show_color_bar:
            font_size_cbar = self._font_size * .6
            cbar = fig.colorbar(im, ax=ax, shrink=0.7)
            cbar.ax.tick_params(labelsize=font_size_cbar)
            if label != "":
                cbar.set_label(label, fontsize=font_size_cbar*1.2)
            else:
                cbar.set_label(self._units, fontsize=font_size_cbar*1.2)

    def plot_slice_eoa(self, c_img, ax, cmap='gray', mode=PlotMode.RASTER, mincbar=np.nan, maxcbar=np.nan) -> None:
        """
        Plots a 2D img for EOA data.
        :param c_img: 2D array
        :param ax: geoaxes
        :return:
        """
        c_ax = ax
        if self._flip_data:
            origin = 'lower'
        else:
            origin = 'upper'

        if self._background == BackgroundType.CARTO_DEF:
            c_ax.stock_img()
        else:
            if self._background == BackgroundType.BLUE_MARBLE_LR:
                img = plt.imread(join(self._eoas_pyutils_path,'viz_utils/imgs/bluemarble.png'))
            if self._background == BackgroundType.BLUE_MARBLE_HR:
                img = plt.imread(join(self._eoas_pyutils_path,'viz_utils/imgs/bluemarble_5400x2700.jpg'))
            if self._background == BackgroundType.TOPO:
                img = plt.imread(join(self._eoas_pyutils_path,'viz_utils/imgs/etopo.png'))
            if self._background == BackgroundType.BATHYMETRY:
                img = plt.imread(join(self._eoas_pyutils_path,'viz_utils/imgs/bathymetry_3600x1800.jpg'))
            c_ax.imshow(img, origin='upper', extent=(-180,180,-90,90), transform=ccrs.PlateCarree())

        if mode == PlotMode.RASTER or mode == PlotMode.MERGED:
            if self._contourf:
                im = c_ax.contourf(self._lons, self._lats, c_img, num_colors=255, cmap='inferno', extent=self._extent)
            else:
                if np.isnan(mincbar):
                    im = c_ax.imshow(c_img, extent=self._extent, origin=origin, cmap=cmap, transform=self._projection)
                else:
                    im = c_ax.imshow(c_img, extent=self._extent, origin=origin, cmap=cmap, vmin=mincbar, vmax=maxcbar, transform=self._projection)

        if mode == PlotMode.CONTOUR or mode == PlotMode.MERGED:
            c_ax.contour(c_img)

        gl = c_ax.gridlines(draw_labels=True, color='grey', alpha=0.5, linestyle='--')
        # gl.xlabel_style = {'size': self._font_size/2, 'color': '#aaaaaa', 'weight':'bold'}
        font_coords = {'size': self._font_size*.6}
        gl.xlabel_style = font_coords
        gl.ylabel_style = font_coords
        gl.top_labels = False
        gl.right_labels = False

        return im

    def get_proper_size(self, rows, cols):
        """
        Obtains the proper size for a figure.
        :param rows: how many rows will the figure have
        :param cols: how many colswill the figure have
        :param prop: Proportion is the proportion to use w/h
        :return:
        """
        if rows == 1:
            return self._figsize * cols * self._fig_prop, self._figsize
        else:
            return self._figsize * cols * self._fig_prop, self._figsize * rows

    def _close_figure(self):
        """Depending on what is disp_images, the figures are displayed or just closed"""
        if self._disp_images:
            plt.show()
        else:
            plt.close()

    def getExtent(self, lats, lons):
        minLat = np.amin(lats)
        maxLat = np.amax(lats)
        minLon = np.amin(lons)
        maxLon = np.amax(lons)
        bbox = (minLon, maxLon, minLat, maxLat)
        return bbox

    def xr_summary(self, ds):
        """ Prints a summary of the netcdf (global attributes, variables, etc)
        :param ds:
        :return:
        """
        print("\n========== Global attributes =========")
        for name in ds.attrs:
            print(F"{name} = {getattr(ds, name)}")

        print("\n========== Dimensions =========")
        for name in ds.dims:
            print(F"{name}: {ds[name].shape}")

        print("\n========== Coordinates =========")
        for name in ds.coords:
            print(F"{name}: {ds[name].shape}")

        print("\n========== Variables =========")
        for cur_variable_name in ds.variables:
            cur_var = ds[cur_variable_name]
            print(F"{cur_variable_name}: {cur_var.dims} {cur_var.shape}")
            
    def nc_summary(self, ds):
        """ Prints a summary of the netcdf (global attributes, variables, etc)
        :param ds: 
        :return: 
        """
        print("\n========== Global attributes =========")
        for name in ds.ncattrs():
            print(F"{name} = {getattr(ds, name)}")

        print("\n========== Variables =========")
        netCDFvars = ds.variables
        for cur_variable_name in netCDFvars.keys():
            cur_var = ds.variables[cur_variable_name]
            print(F"Dimensions for {cur_variable_name}: {cur_var.dimensions} {cur_var.shape}")

    def add_roads(self, ax):
        # Names come from: https://www.naturalearthdata.com/features/
        # -- Add states
        roads = cfeature.NaturalEarthFeature(
            category='cultural',
            name='roads',
            scale='10m',
            facecolor='none')

        ax.add_feature(roads, edgecolor='black')
        return ax

    def add_states(self, ax):
        # Names come from: https://www.naturalearthdata.com/features/
        # -- Add states
        states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')

        ax.add_feature(states_provinces, edgecolor='gray')
        return ax

    def plot_3d_data_npdict(self, np_variables:list, var_names:list, z_levels= [], title='',
                          file_name_prefix='', cmap='viridis', z_names = [], flip_data=True,
                          show_color_bar=True, plot_mode=PlotMode.RASTER, mincbar=np.nan, maxcbar=np.nan):
        """
        Plots multiple z_levels for multiple fields.
        It uses rows for each depth, and columns for each variable
        """
        create_folder(self._output_folder)

        # If the user do not requires any z-leve, then all are plotted
        if len(z_levels) == 0:
            z_levels = range(np_variables[0].shape[0])

        rows = len(z_levels)
        cols = len(var_names)

        fig, _axs = plt.subplots(rows, cols,
                                 figsize=self.get_proper_size(rows, cols),
                                 subplot_kw={'projection': self._projection})

        for c_row, c_slice in enumerate(z_levels):  # Iterates over the z-levels
            # Verify the index of the z_levels are the original ones.
            if len(z_names) != 0:
                c_slice_txt = z_names[c_slice]
            else:
                c_slice_txt = c_slice

            c_mincbar = np.nan
            c_maxcbar = np.nan
            for idx_var, c_var in enumerate(var_names): # Iterate over the fields
                if rows*cols == 1:
                    ax = _axs
                else:
                    ax = _axs.flatten()[c_row*len(var_names)+ idx_var]

                # Here we chose the min and max colorbars for each field
                if not(np.all(np.isnan(mincbar))):
                    c_mincbar = mincbar[idx_var]
                if not(np.all(np.isnan(maxcbar))):
                    c_maxcbar = maxcbar[idx_var]

                if self._auto_colormap:
                    cmap = select_colormap(c_var)

                im = self.plot_slice_eoa(np_variables[c_var][c_slice,:,:], ax, cmap=cmap, mode=plot_mode,
                                         mincbar=c_mincbar, maxcbar=c_maxcbar)

                if self._show_var_names:
                    c_title = F'{var_names[idx_var]} {title} Z-level:{c_slice_txt}'
                else:
                    c_title = F'{title} Z-level:{c_slice_txt}'

                ax.set_title(c_title, fontsize=self._font_size)

                self.add_colorbar(fig, im, ax, show_color_bar)

        file_name = F'{file_name_prefix}_{c_slice_txt:04d}'
        pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
        self._close_figure()

    def plot_2d_data_np(self, np_variables:list, var_names:list, title='',
                            file_name_prefix='', cmap='viridis',  flip_data=True,
                            rot_90=True, show_color_bar=True, plot_mode=PlotMode.RASTER, mincbar=np.nan, maxcbar=np.nan):
        '''
        Wrapper function to receive raw 2D numpy data. It calls the 'main' function for 3D plotting
        :param np_variables:
        :param var_names:
        :param title:
        :param file_name_prefix:
        :param cmap:
        :param flip_data:
        :param rot_90:
        :param show_color_bar:
        :param plot_mode:
        :param mincbar:
        :param maxcbar:
        :return:
        '''
        npdict_3d = {field_name: np.expand_dims(np_variables[i,:,:]) for i, field_name in enumerate(var_names)}
        self.plot_3d_data_npdict(npdict_3d, var_names, z_levels=[0], title=title,
                            file_name_prefix=file_name_prefix, cmap=cmap, z_names = [], flip_data=flip_data,
                            show_color_bar=show_color_bar, plot_mode=plot_mode, mincbar=mincbar, maxcbar=maxcbar)

