# %%
from download_data.Download_SSS_SMAP_satellite import parallel_sss_download
from multiprocessing import Pool

TOT_PROC = 1
# %%

# %% ----------- Download SSS data
output_folder = "/Net/work/ozavala/DATA/GOFFISH/SSS/SMAP_Global"
# output_folder = "../Data/AVISO/SSS/SMAP_Global"
years = range(2013,2022)

parallel_sss_download(output_folder, years)

# %% ----------- Download SSS data parallel
# p = Pool(TOT_PROC)
# p.map(parallel_sss_download, range(TOT_PROC))