import numpy as np
import xarray as xr
import os
import sys

# read in data and equi files
fname = 'ARGO_WO_2024_PFL_D_labsea_splitcost'
fnames = [fname, 'ad'+fname]

iter = sys.argv[1]
ext = sys.argv[2]

root_dir = '/scratch/08382/shoshi/labsea_runs/'

run_dir = root_dir + 'assim_argo_MG/run_adhi_it' + iter + ext + '/'
ds_argo = xr.open_dataset(run_dir + fname + '.nc')

hr_files = sorted(glob.glob(adj_dir + 'profiles/' + fname + '.*.equi.nc'))

def prepare_hr(ds):
    # Ensure iPROF uses the unique global index so concatenation works
    return ds.assign_coords(iPROF=ds.prof_ind_glob.values)

# Load all equi files into one dataset
print(f"Loading {len(hr_files)} files into hr_all...")
hr_all = xr.open_mfdataset(
    hr_files, 
    concat_dim="iPROF", 
    combine="nested", 
    preprocess=prepare_hr,
    chunks={'iPROF': 100} # Optional: helps with memory if dataset is massive
)

# Subset data to the active profiles in runtime
hr_glob_indices = np.unique(hr_all.prof_ind_glob.values.astype(int))
argo_sub = ds_argo.isel(iPROF=hr_glob_indices - 1).copy()

# Extract Lo-Res Grid Arrays
grid_dir_hr = root_dir + 'grid_hires/'
grid_dir_lr = root_dir + 'grid_lores/'

grid_hr = xmitgcm.open_mdsdataset(grid_dir_hr, grid_dir=grid_dir_hr, iters = None)
grid_lr = xmitgcm.open_mdsdataset(grid_dir_lr, grid_dir=grid_dir_lr, iters = None)

hfacc_hr = rdmds(root_dir + 'grid_hires_cleanbathy/' + 'hFacC')
hfacc_lr = rdmds(root_dir + 'grid_lores_cleanbathy/'+ 'hFacC')

grid_hr.hFacC.values = hfacc_hr
grid_lr.hFacC.values = hfacc_lr

xc_lr = grid_lr.XC.values
yc_lr = grid_lr.YC.values
tile_size = 24

# For the "End" of the tile (sNx+1 logic), we need the value 
# one index beyond the tile. If it's the last tile, we estimate the next center.
dx = 0.3333333333
dy = 0.3333333333

lons = argo_sub.prof_lon.values
lats = argo_sub.prof_lat.values

assigned_tile_id = np.zeros(len(argo_sub.iPROF), dtype=int)

for j in range(npy_lr):
    for i in range(npx_lr):
        tile_id = j * npx_lr + i + 1
        
        # Start coordinate (1,1)
        x1 = xc_lr[i * tile_size]
        y1 = yc_lr[j * tile_size]
        
        # End coordinate logic (sNx+1)
        # If at the very edge of the global grid, we extrapolate the next center
        if (i + 1) * tile_size < len(xc_lr):
            x_end = xc_lr[(i + 1) * tile_size]
        else:
            x_end = xc_lr[-1] + dx
            
        if (j + 1) * tile_size < len(yc_lr):
            y_end = yc_lr[(j + 1) * tile_size]
        else:
            y_end = yc_lr[-1] + dy

        # Replicate the Fortran check: [x1 <= Lon < x_end] AND [y1 <= Lat < y_end]
        # Includes the 360 wrap-around logic from your snippet
        mask_standard = (lons >= x1) & (lons < x_end) & (lats >= y1) & (lats < y_end)
        
        # Wrap-around check (if grid crosses the date line/zero)
        mask_wrap = (x_end < x1) & (lons + 360 >= x1) & (lons + 360 < x_end + 360) & \
                    (lats >= y1) & (lats < y_end)
        
        assigned_tile_id[mask_standard | mask_wrap] = tile_id

# Map to Ocean-Only and Attach
prof_lr_ocean_num = np.array([lr_ocean_map.get(t, -999) for t in assigned_tile_id])
argo_sub['prof_LR_tile'] = (['iPROF'], prof_lr_ocean_num)
argo_sub['prof_ind_glob'] = (['iPROF'], hr_glob_indices)

# 4. Final Data Mapping (Equi data)
hr_all_indexed = hr_all.set_index(iPROF='prof_ind_glob')
argo_sub['prof_Tequi'] = (['iPROF', 'iDEPTH'], hr_all_indexed.prof_T.sel(iPROF=argo_sub.prof_ind_glob).values)
argo_sub['prof_Sequi'] = (['iPROF', 'iDEPTH'], hr_all_indexed.prof_S.sel(iPROF=argo_sub.prof_ind_glob).values)

argo_sub['prof_Tmask'] = (['iPROF', 'iDEPTH'], hr_all_indexed.prof_Tmask.sel(iPROF=argo_sub.prof_ind_glob).values)
argo_sub['prof_Smask'] = (['iPROF', 'iDEPTH'], hr_all_indexed.prof_Smask.sel(iPROF=argo_sub.prof_ind_glob).values)


# Define the output directory
out_dir = os.path.join(run_dir, 'prof_LR_equi')
os.makedirs(out_dir, exist_ok=True)

for f_base in fnames:
    # Group the subsetted dataset by the Lo-Res ocean tile number
    grouped = argo_sub.groupby('prof_LR_tile')

    for tile_num, tile_ds in grouped:
        if tile_num <= 0:  # Skip land (-999) or unassigned (0)
            continue

        # Drop the original observation variables to avoid the name conflict
        # We also drop other non-essential variables to keep the equi file small
        tmp_ds = tile_ds.drop_vars(['prof_T', 'prof_S'], errors='ignore')

        # Rename the equivalent variables to the standard MITprof names
        export_ds = tmp_ds.rename({
            'prof_Tequi': 'prof_T',
            'prof_Sequi': 'prof_S'
        })

        # Select only the variables required for the MITgcm equi file
        # This keeps the file structure identical to what MITgcm produces
        vars_to_keep = ['prof_ind_glob', 'prof_T', 'prof_S', 'prof_Tmask', 'prof_Smask']
        export_ds = export_ds[vars_to_keep]

        # Sort by profile index and set output path
        export_ds = export_ds.sortby('prof_ind_glob')
        out_file = os.path.join(out_dir, f"{fname}.{int(tile_num):03d}.001.equi.nc")

        # Write to NetCDF
        export_ds.to_netcdf(out_file, mode='w', format='NETCDF4_CLASSIC')

print(f"Successfully wrote Lo-Res equi files to: {out_dir}")
