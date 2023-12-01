# %%
from download_data.Download_SSS_SMAP_satellite import parallel_sss_download
from multiprocessing import Pool

TOT_PROC = 10
# %%

# %% ----------- Download SSS data
output_folder = "/Net/work/ozavala/DATA/GOFFISH/SSS/SMAP_Global"
# output_folder = "../Data/AVISO/SSS/SMAP_Global"
years = range(2013,2016)
# For testing single proc
# parallel_sss_download(output_folder, years)
p = Pool(TOT_PROC)
p.starmap(parallel_sss_download, [(output_folder, years, x, TOT_PROC) for x in range(TOT_PROC)])