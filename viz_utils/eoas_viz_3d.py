from os.path import join

import SimpleITK as sitk

import plotly
import plotly.figure_factory as ff
import plotly.graph_objs as go
import numpy as np
from skimage import measure

class ImageVisualizer3D:
    _output_folder = 'output_medical'
    _open_browser = False
    _COLORS = ['y', 'r', 'c', 'b', 'g', 'w', 'k', 'y', 'r', 'c', 'b', 'g', 'w', 'k']

    def __init__(self, **kwargs):
        # All the arguments that are passed to the constructor of the class MUST have its name on it.
        for arg_name, arg_value in kwargs.items():
            self.__dict__["_" + arg_name] = arg_value

    def __getattr__(self, attr):
        '''Generic getter for all the properties of the class'''
        return self.__dict__["_" + attr]

    def __setattr__(self, attr, value):
        '''Generic setter for all the properties of the class'''
        self.__dict__["_" + attr] = value

    def _check_file_name(self, file_name):
        """It only validates that the name contains html at the end"""
        return file_name if file_name.find('html') != -1 else F'{file_name}.html'

    def plot_surface_itk(self, ctr_np, title='', file_name=''):
        """
        Simple wrapper for itk image formats
        :param ctr_np:
        :param title:
        :param file_name:
        :return:
        """
        self.plot_surface_np(sitk.GetArrayFromImage(ctr_np), title=title, file_name=file_name)


    def plot_surface_np(self, ctr_np, title='', file_name=''):
        """
        Makes a 3D surface using Ploty from a 3D Volume. It expects a binary mask
        :param ctr_np:
        :param title:
        :param file_name:
        :return:
        """
        # Computes marching cubes to obtain a mesh
        print('\tMarching cubes...')
        vertices, simplices, normals, values = measure.marching_cubes_lewiner(ctr_np,.8)
        print('\tDone!')

        print('\t3D Surface...')
        x, y, z = zip(*vertices)
        camera = dict(up=dict(x=1, y=0, z=0),
                      center=dict(x=0, y=0, z=0),
                      eye=dict(x=2, y=2, z=0.1))

        fig = ff.create_trisurf(x=y,
                                y=z,
                                z=x,
                                plot_edges=True,
                                color_func=None,
                                gridcolor='rgb(0,0,0)',
                                simplices=simplices,
                                show_colorbar=False,
                                title=title)

        fig['layout'].update(
            scene=dict(camera=camera),
            title=title)

        plotly.offline.plot(fig, filename=join(self.output_folder,self._check_file_name(file_name)), auto_open=self.open_browser)


    def plot_scatter_np(self, ctr_np, file_name, title=''):
        """
        Makes a 3D scatter plot with the obtained
        :param ctr_np:
        :param file_name:
        :param title:
        :return:
        """
        z, y, x = np.where(ctr_np > 0)
        data=[]
        max_z = max(z)
        min_z= min(z)
        min_y= min(y)
        min_x= min(x)

        # Making data and assigning colors
        print("\tMaking data")
        for z_level in range(int(min_z), int(max_z)):
            color = 'rgb({},0,0)'.format(int(255*(z_level-min_z)/(max_z-min_z)))
            idx = np.where(z == z_level)

            data.append(go.Scatter3d(x=x[idx]-min_x, y=y[idx]-min_y, z=z[idx]-min_z, hoverinfo=None,
                                     mode='markers', marker=dict(symbol='circle',size=2, color=color)))

        # =========== ALL THIS PART IS TO 'IMPROVE' THE LAYOUT ===========
        tick_l = 5
        tick_w = 3
        gridwidth = tick_w
        zerolinewidth = tick_w + 3
        zero_color = 'rgb(0,0,0)'
        grid_color = 'rgb(150,150,150)'
        tickfont = dict( color='black', size=19, family='Old Standard TT, serif', )
        titlefont = dict( color='black', size=35, family='Old Standard TT, serif', )

        # https://plot.ly/python/3d-axes/:w
        layout = go.Layout(
            scene=dict(
                xaxis=dict(
                    nticks=10,
                    gridwidth=gridwidth,
                    tickfont=tickfont,
                    showbackground=True,
                    # zerolinecolor=zero_color,
                    # zerolinewidth=zerolinewidth,
                    backgroundcolor='rgb(255,255,255)',
                    # zeroline=False,
                    gridcolor=grid_color,
                    # ticklen=tick_l,
                    # tickwidth=tick_w,
                    titlefont=titlefont
                ),
                yaxis=dict(
                    nticks=5,
                    gridwidth=gridwidth,
                    showbackground=True,
                    tickfont=tickfont,
                    # zeroline=False,
                    # zerolinewidth=zerolinewidth,
                    # zerolinecolor=zero_color,
                    backgroundcolor='rgb(255,255,255)',
                    gridcolor=grid_color,
                    # ticklen=tick_l,
                    # tickwidth=tick_w,
                    titlefont=titlefont
                ),
                zaxis=dict(
                    nticks=5,
                    gridwidth=gridwidth,
                    showbackground=True,
                    tickfont=tickfont,
                    # zeroline=False,
                    # zerolinewidth=zerolinewidth,
                    # zerolinecolor=zero_color,
                    backgroundcolor='rgb(255,255,255)',
                    gridcolor=grid_color,
                    # ticklen=tick_l,
                    # tickwidth=tick_w,
                    titlefont=titlefont),
                aspectmode='cube',
                camera=dict(
                    up=dict(x=1, y=0, z=0),
                    center=dict(x=0, y=0, z=0),
                    eye=dict(x=2, y=-3, z=1.5))
            ))

        fig = go.Figure(data=data, layout=layout)

        file_name = file_name if file_name.find('html') != -1 else F'{file_name}.html'
        fig['layout'].update( title=title)
        plotly.offline.plot(fig, filename=join(self.output_folder,self._check_file_name(file_name)), auto_open=self.open_browser)
