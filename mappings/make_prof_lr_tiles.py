import numpy as np
import xarray as xr
import xmitgcm
import glob
import os

import sys
sys.path.append('/home/shoshi/MITgcm_c68r/MITgcm/utils/python/MITgcmutils')
from MITgcmutils import rdmds

sys.path.append('/home/shoshi/jupyter_notebooks')
from read_write import *

sys.path.append('/home/shoshi/jupyter_notebooks/DivA/swot_assim/')
import map_ocean_tiles

fname = 'ARGO_WO_2024_PFL_D_labsea_splitcost'
fnames = [fname, 'ad'+fname]

iter = sys.argv[1]
ext = sys.argv[2]

run_dir = '/scratch/shoshi/labsea_MG_12/assim_argo_MG/run_adhi_it' + iter + ext + '/'
run_dir_lo = '/scratch/shoshi/labsea_MG_12/assim_argo_MG/run_adlo_it' + iter + ext + '/'

grid_dir_hr = '/scratch/shoshi/labsea_MG_12/grid_hires/'
grid_dir_lr = '/scratch/shoshi/labsea_MG_12/grid_lores/'

grid_hr = xmitgcm.open_mdsdataset(grid_dir_hr, grid_dir=grid_dir_hr, iters = None)
grid_lr = xmitgcm.open_mdsdataset(grid_dir_lr, grid_dir=grid_dir_lr, iters = None)

hfacc_hr = rdmds('/scratch/shoshi/labsea_MG_12/grid_hires_cleanbathy/' + 'hFacC')
#hfacc_lr = rdmds(run_dir_lo + 'hFacC')
hfacc_lr = rdmds('/scratch/shoshi/labsea_MG_12/grid_lores_cleanbathy/'+ 'hFacC')

grid_hr.hFacC.values = hfacc_hr
grid_lr.hFacC.values = hfacc_lr

npx_lr = 5
npy_lr = 4

npx_hr = 20
npy_hr = 16

# create mapping of all tiles to ocean tiles
lr_mapping = map_ocean_tiles.tile_to_ocean_number(grid_lr.hFacC[0].values, npx_lr, npy_lr)
hr_mapping = map_ocean_tiles.tile_to_ocean_number(grid_hr.hFacC[0].values, npx_hr, npy_hr)

# create mapping of high-res ocean tiles in low-res ocean tiles
hi_to_lo_mapping = map_ocean_tiles.map_hr_to_lr_by_ocean_number(lr_mapping, hr_mapping, npx_lr, npy_lr, npx_hr, npy_hr)

def get_hr_ocean_files(run_dir, f, ocean_mapping, lr_ocean_tile):
    """
    Return a list of HR ocean files corresponding to a given LR ocean tile.
    """
    hr_tiles = ocean_mapping.get(lr_ocean_tile, [])
    file_list = []
    for tile in hr_tiles:
        file_path = os.path.join(run_dir + 'profiles/', f"{f}.{tile:03d}.001.equi.nc")
        if os.path.exists(file_path):
            file_list.append(file_path)
    return file_list

def combine_equi_files(files):

    datasets = []
    for f in files:
        ds = xr.open_dataset(f)
        datasets.append(ds)

    # Concatenate along iPROF (even if lengths differ)
    equi = xr.concat(datasets, dim='iPROF', combine_attrs='override')

    # Rebuild a continuous iPROF coordinate
    equi = equi.assign_coords(iPROF=range(equi.dims['iPROF']))
    equi = equi.sortby('prof_ind_glob')

    return equi


out_dir = os.path.join(run_dir, 'prof_LR_equi')
os.makedirs(out_dir, exist_ok=True)

for f in fnames:
    for tile in hi_to_lo_mapping.keys():

        files_tile = get_hr_ocean_files(run_dir , f, hi_to_lo_mapping, lr_ocean_tile=tile)
        if files_tile :
            tile_combined = combine_equi_files(files_tile)
            out_file = os.path.join(out_dir, f"{f}.{tile:03d}.001.equi.nc")
            tile_combined.to_netcdf(out_file)
