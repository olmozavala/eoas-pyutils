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

def read_multiple_netcdf_xarr(file_names:str, fields: list):
    ds = []
    for c_file in file_names:
        ds = xr.load_dataset(c_file)

