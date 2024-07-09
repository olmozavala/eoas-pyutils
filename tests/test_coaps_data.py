# Testing the provided functions
import pytest
import matplotlib.pyplot as plt
from io_utils.coaps_io_data import get_aviso_by_date, get_aviso_by_month, get_sss_by_date, get_sst_ghrsst_by_date, get_chlora_noaa_by_date
from io_utils.dates_utils import get_day_of_year_from_month_and_day 
from datetime import datetime

c_date = datetime(2016, 1, 10)
c_date_str = c_date.strftime("%Y-%m-%d")
bbox = [17.5, 32.5, -98, -76]

def test_aviso():
    aviso_folder = "/unity/f1/ozavala/DATA/GOFFISH/AVISO/GoM/"
    # Print 
    print(f"Reading monthly aviso data... ")
    aviso_data, lats, lons = get_aviso_by_month(aviso_folder, c_date, bbox)
    print(f"The number of times available are: {aviso_data.time.size} for date {c_date_str}")
    assert aviso_data.time.size == 31
    assert lats.size == 60
    assert aviso_data.adt.shape == (31, 60,88)

    print(f"Reading daily aviso data... ")
    aviso_data, lats, lons = get_aviso_by_date(aviso_folder, c_date, bbox)
    print(f"The number of times available are: {aviso_data.time.size} requested {c_date_str} available {aviso_data.time.values}")
    assert aviso_data.time.size == 1
    assert lats.size == 60
    assert aviso_data.adt.shape == (60,88)
    # plt.imshow(aviso_data.adt[:,:], origin="lower")
    # plt.show()

def test_satellite_sst():
    sst_folder = "/unity/f1/ozavala/DATA/GOFFISH/SST/OISST"

    print(f"Reading SST data... ")
    sst_data, lats, lons = get_sst_ghrsst_by_date(sst_folder, c_date, bbox)
    assert sst_data.analysed_sst.shape == (1, 1351, 2201)
    assert lats.size == (1351)
    assert lons.size == (2201)
    # plt.imshow(sst_data.analysed_sst[0,:,:], origin="lower")
    # plt.show()

def test_satellite_sss():
    sss_folder = "/unity/f1/ozavala/DATA/GOFFISH/SSS/SMAP_Global"
    c_date = datetime(2016, 1, 10)
    c_date_str = c_date.strftime("%Y-%m-%d")

    print(f"Reading SSS data... ")
    sss_data, lats, lons = get_sss_by_date(sss_folder, c_date, bbox)
    print(sss_data.sss_smap.shape)
    assert sss_data.sss_smap.shape == (60,88)
    assert lats.size == (60)
    assert lons.size == (88)
    # plt.imshow(sss_data.sss_smap[:,:], origin="lower")
    # plt.show()

def test_satellite_chlora_noaa():
    chlora_folder = "/unity/f1/ozavala/DATA/GOFFISH/CHLORA/NOAA"
    c_date = datetime(2019, 1, 10)

    print(f"Reading chlora data... ")
    chlora_data, lats, lons = get_chlora_noaa_by_date(chlora_folder, c_date, bbox)
    print(chlora_data.variables)
    assert chlora_data.chlor_a.squeeze().shape == (180, 264)
    assert lats.size == (180)
    assert lons.size == (264)
    # plt.imshow(chlora_data.chlor_a.squeeze())
    # plt.show()