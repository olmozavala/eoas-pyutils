# %%
import os
import shutil
# from download_data.Download_SSS_SMAP_satellite import parallel_sss_download
from ..download_data.Download_SSS_SMAP_satellite import parallel_sss_download  # If running from the root folder

# %%
def test_parallel_sss_download():
    # Setup
    proc_id = 1
    output_folder = "test_folder"
    years = [2020]

    # Execute
    parallel_sss_download(proc_id, output_folder, years)

    # Verify
    for year in years:
        year_folder = os.path.join(output_folder, str(year))
        assert os.path.exists(year_folder), f"Year folder {year_folder} does not exist"
        for day in range(1, 366):
            if day % TOT_PROC == proc_id:
                file_name = f"RSS_smap_SSS_L3_8day_running_{year}_{day:03d}_FNL_v05.0.nc"
                file_path = os.path.join(year_folder, file_name)
                assert os.path.exists(file_path), f"File {file_path} does not exist"

    # Cleanup
    shutil.rmtree(output_folder)

test_parallel_sss_download()