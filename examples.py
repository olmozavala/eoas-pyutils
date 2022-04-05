import sys
sys.path.append("hycom_utils/python")

from hycom.io import read_hycom_fields, subset_hycom_field, read_hycom_coords
from hycom.info import read_field_names
from viz_utils.eoa_viz import EOAImageVisualizer

## Reading data
print("Reading data...")
input_file = "./hycom_utils/test_data/archv.2009_153_00.a"
coords_file = "./hycom_utils/test_data/regional.grid.a"
layers = [0, 1, 2, 3]  # Depth layers we want to read
print(F"The fields available are: {read_field_names(input_file)}")
# Reading specific field and layers
fields = ['srfhgt', 'temp', 'u-vel.']   # Which fields to plot
hycom_fields = read_hycom_fields(input_file, fields, layers)

print(F"The coords available are: {read_field_names(coords_file)}")
coords = read_hycom_coords(coords_file, ['plon:', 'plat:'])
lons = coords['plon:']
lats = coords['plat:']
print("Done!")

## Plotting data
viz_obj = EOAImageVisualizer(lats=lats, lons=lons, disp_images=True, output_folder="outputs", eoas_pyutils_path=".")
viz_obj.plot_3d_data_npdict(hycom_fields, ['temp','u-vel.'], range(len(layers)), 'MyTitle', 'myplot')
print("Done!")
##

