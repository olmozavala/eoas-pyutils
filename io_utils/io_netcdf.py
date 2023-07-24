from enum import Enum
import os
from os import walk, listdir
from os.path import join

import numpy as np
from netCDF4 import Dataset
import xarray as xr

def read_netcdf_xr(file_name:str, fields: list,  replace_to_nan=True, rename_fields=[]):
    nc_file = xr.load_dataset(file_name)
    all_fields =  list(nc_file.variables)

    if len(fields) == 0:
        fields = all_fields
        # print(F"Reading all the fields in the file: {fields}")
    if not(np.all([field in all_fields for field in fields])):
        print(F"Warning!!!!! Fields {[field for field in fields if not(field in all_fields)]} are not"
              F" in the netcdf file {file_name}, removing them from the list.")
        fields = [field for field in fields if field in all_fields ]

    # This is just a patch to 'rename' the variables on the fly
    if len(rename_fields) > 0:
        nc_fields = {rename_fields[idx]: all_fields[field] for idx, field in enumerate(fields)}
    else:
        nc_fields = {field: nc_file[field] for field in fields}

    return nc_fields


def read_netcdf(file_name:str, fields: list,  replace_to_nan=True, rename_fields=[]):
    nc_file = Dataset(file_name, "r", format="NETCDF4")
    all_fields =  nc_file.variables

    if len(fields) == 0:
        fields = all_fields
        # print(F"Reading all the fields in the file: {fields}")
    if not(np.all([field in all_fields for field in fields])):
        print(F"Warning!!!!! Fields {[field for field in fields if not(field in all_fields)]} are not"
              F" in the netcdf file {file_name}, removing them from the list.")
        fields = [field for field in fields if field in all_fields ]

    # This is just a patch to 'rename' the variables on the fly
    if len(rename_fields) > 0:
        nc_fields = {rename_fields[idx]: all_fields[field] for idx, field in enumerate(fields)}
    else:
        nc_fields = {field: all_fields[field] for field in fields}

    return nc_fields

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

def read_multiple_netcdf_xarr(file_names:str, fields: list):
    ds = []
    for c_file in file_names:
        ds = xr.load_dataset(c_file)


# This is an enum class that provides the standar names used for differnet variables
class CF_StandardNames(Enum):
    # Standard names for the variables
    # http://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
    TEMP = "sea_water_temperature"
    SALINITY = "sea_water_salinity"
    SSH = "sea_surface_height"
    U = "eastward_sea_water_velocity"
    V = "northward_sea_water_velocity"
    MLD = "sea_water_mixed_layer_thickness"
    CHLORA = "chlorophyll_concentration_in_sea_water"
