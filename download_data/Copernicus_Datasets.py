

# Make an enumeration of the datasets
from enum import Enum
import datetime

class Copernicus_Enum(Enum):
    CHLORA_L3_D = 1 # https://data.marine.copernicus.eu/product/OCEANCOLOUR_GLO_BGC_L3_MY_009_103/description
    SST_ODYSSEA_L3_D = 2 # https://data.marine.copernicus.eu/product/SST_GLO_SST_L3S_NRT_OBSERVATIONS_010_010/description
    SST_OSTIA_L4_D = 3 # https://data.marine.copernicus.eu/product/SST_GLO_SST_L4_NRT_OBSERVATIONS_010_001/description

# Create a dictionary with the Copernicus datasets
Copernicus_Datasets = {
    Copernicus_Enum.CHLORA_L3_D: {
        'name': 'Global Ocean Color Plankton and Reflectances MY L3 daily observations',
        'short_name': 'Ocean_Color',
        'id': 'c3s_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D',
        'version': '202303',
        'start_date': '?',
        'resolution': '4km 0.036',
        'variables': ["CHL", "MICRO", "MICRO_BIAS", "MICRO_RMSE", "NANO", "NANO_BIAS", "NANO_RMSE", "PICO", "PICO_BIAS", "PICO_RMSE"],
    },
    Copernicus_Enum.SST_ODYSSEA_L3_D: {
        'name': 'ODYSSEA Global Ocean - Sea Surface Temperature Multi-sensor L3 Observations',
        'short_name': 'ODYSSEA_SST',
        'id': 'IFREMER-GLOB-SST-L3-NRT-OBS_FULL_TIME_SERIE',
        'version': '202211',
        'start_date': '2021-01-01',
        'resolution': '0.1x0.1',
        'variables': ["adjusted_sea_surface_temperature", "bias_to_reference_sst", "or_latitude", "or_longitude", "sea_surface_temperature", "sst_dtime"],
    },
    Copernicus_Enum.SST_OSTIA_L4_D: {
        'name': 'Global Ocean OSTIA Sea Surface Temperature and Sea Ice Analysis',
        'short_name': 'OSTIA_SST',
        'id': 'METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2',
        'version': '',
        'start_date': '2007-01-01',
        'resolution': '0.05x0.05',
        'variables': ["analysed_sst", "analysis_error", "mask", "sea_ice_fraction"],
    },
}
