import numpy as np
import xarray as xr
import os
import sys
import glob
import xmitgcm

sys.path.append('/work2/08382/shoshi/stampede3/MITgcm_c69j/MITgcm/utils/python/MITgcmutils')
from MITgcmutils import rdmds

fname = 'swot_obsfit_cycles_9thru11_labsea_L3v3'
fnames = [fname, 'ad'+fname]

iter = sys.argv[1]
ext = sys.argv[2]
rundirs = sys.argv[3]

root_dir = '/scratch/08382/shoshi/labsea_runs/'

run_dir = rundirs + '/run_adhi_it' + iter + ext + '/'

forward_files = sorted(glob.glob(f"{run_dir}{fname}.*.equi.nc"))
all_tile_datasets = []

### combine all equi files
for f_path in forward_files:
    # Extract the tile string (e.g., '003.001') from the filename
    tile_suffix = ".".join(os.path.basename(f_path).split('.')[1:3])
    ad_path = os.path.join(run_dir , f"ad{fname}.{tile_suffix}.equi.nc")

    # Load Forward
    ds_fwd = xr.open_dataset(f_path)

    if os.path.exists(ad_path):
        # Load Adjoint and rename variables
        ds_ad = xr.open_dataset(ad_path)
        rename_dict = {v: f"ad{v}" for v in ds_ad.data_vars if v != 'sample_ind_glob'}
        ds_ad = ds_ad.rename(rename_dict)

        # Merge on iSAMPLE (since they share local indexing within the same high-res tile)
        # We ensure they align on sample_ind_glob 
        combined_tile = xr.merge([ds_fwd, ds_ad.drop_vars('sample_ind_glob', errors='ignore')])
    else:
        # If adjoint is missing for this tile, fill with zeros to maintain structure
        combined_tile = ds_fwd.copy()
        for v in ds_fwd.data_vars:
            if v != 'sample_ind_glob':
                combined_tile[f"ad{v}"] = xr.zeros_like(ds_fwd[v])

    # Use the global index as the coordinate for the upcoming concatenation
    combined_tile = combined_tile.assign_coords(iSAMPLE=combined_tile.sample_ind_glob.values.astype(int))
    all_tile_datasets.append(combined_tile)

# Combine everything into one master dataset
hr_all = xr.concat(all_tile_datasets, dim='iSAMPLE').sortby('iSAMPLE')

### load grid info
grid_dir_lr = root_dir + 'grid_lores/'
grid_lr = xmitgcm.open_mdsdataset(grid_dir_lr, grid_dir=grid_dir_lr, iters=None)
hfacc_lr = rdmds(root_dir + 'grid_lores_cleanbathy/' + 'hFacC')
grid_lr.hFacC.values = hfacc_lr

xc_lr = rdmds(grid_dir_lr + 'XC')
yc_lr = rdmds(grid_dir_lr + 'YC')
ny, nx = xc_lr.shape
tile_size = 24

npx_lr = 5
npy_lr = 4

npx_hr = 20
npy_hr = 16


### Build lr_ocean_map
surface_mask = hfacc_lr[0, :, :]
lr_ocean_map = {}
ocean_count = 0
for j in range(npy_lr):
    for i in range(npx_lr):
        tile_id = j * npx_lr + i + 1
        if np.any(surface_mask[j*tile_size:(j+1)*tile_size, i*tile_size:(i+1)*tile_size] > 0):
            ocean_count += 1
            lr_ocean_map[tile_id] = ocean_count


### Map to Lo-Res Tiles
ds_argo = xr.open_dataset(run_dir + fname + '.nc')
hr_indices = hr_all.iSAMPLE.values
argo_sub = ds_argo.isel(iSAMPLE=hr_indices - 1).copy()

lons, lats = argo_sub.sample_lon.values, argo_sub.sample_lat.values
dx = dy = 0.3333333333
assigned_tile_id = np.zeros(len(argo_sub.iSAMPLE), dtype=int)

for j in range(npy_lr):
    for i in range(npx_lr):
        tile_id = j * npx_lr + i + 1
        x1, y1 = xc_lr[j*tile_size, i*tile_size], yc_lr[j*tile_size, i*tile_size]
        x_end = xc_lr[j*tile_size, (i+1)*tile_size] if (i+1)*tile_size < nx else xc_lr[0, -1] + dx
        y_end = yc_lr[(j+1)*tile_size, i*tile_size] if (j+1)*tile_size < ny else yc_lr[-1, 0] + dy

        mask = (lons >= x1) & (lons < x_end) & (lats >= y1) & (lats < y_end)
        wrap = (x_end < x1) & (lons + 360 >= x1) & (lons + 360 < x_end + 360) & (lats >= y1) & (lats < y_end)
        assigned_tile_id[mask | wrap] = tile_id

# Attach metadata to hr_all for grouping
hr_all['obs_LR_tile'] = (['iSAMPLE'], np.array([lr_ocean_map.get(t, -999) for t in assigned_tile_id]))

### Group and Write Out
out_dir = os.path.join(run_dir, 'obs_LR_equi')
os.makedirs(out_dir, exist_ok=True)

vars_to_keep = ['sample_ind_glob', 'mod_val', 'mod_mask']
grouped = hr_all.groupby('obs_LR_tile')

for tile_num, tile_ds in grouped:
    if tile_num <= 0: continue

    tile_str = f"{int(tile_num):03d}.001"

    # --- Write Forward File ---
    fwd_out = tile_ds[vars_to_keep].sortby('sample_ind_glob')
    fwd_out.to_netcdf(os.path.join(out_dir, f"{fname}.{tile_str}.equi.nc"), format='NETCDF4_CLASSIC')

    # --- Write Adjoint File ---
    ad_vars_to_extract = [f"ad{v}" for v in vars_to_keep if v != 'sample_ind_glob']
    rename_map = {f"ad{v}": v for v in vars_to_keep if v != 'sample_ind_glob'}

    ad_out = tile_ds[['sample_ind_glob'] + ad_vars_to_extract].rename(rename_map)
    ad_out = ad_out.sortby('sample_ind_glob')

    ad_out.to_netcdf(os.path.join(out_dir, f"ad{fname}.{tile_str}.equi.nc"), format='NETCDF4_CLASSIC')


