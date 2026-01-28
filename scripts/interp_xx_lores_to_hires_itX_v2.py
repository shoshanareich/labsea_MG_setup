import numpy as np
from scipy.interpolate import griddata, Rbf 
import sys
sys.path.append('/home/shoshi/MITgcm_c68r/MITgcm/utils/python/MITgcmutils')
from MITgcmutils import rdmds

nx=120
ny=96
nz=21

factor = 4 # lowres * factor = hires

nxh=nx*factor
nyh=ny*factor

#dirroot='/home/shoshi/MITgcm_c68r/MITgcm/verification/lab_sea/'
dirroot='/scratch/shoshi/labsea_MG_12/'

#iter = '0000'
iter = sys.argv[1] # sys.argv[0] is name of python file
print(iter)

griddir_lr = '/scratch/shoshi/labsea_MG_12/grid_lores/'
griddir_hr = '/scratch/shoshi/labsea_MG_12/grid_hires/'

#dirrun_lr = dirroot + 'run_adlo_it' + iter + '/'
#dirrun_hr = dirroot + 'run_adhi_it' + iter + '/'

dirrun_pup = dirroot + 'assim_argo_MG/run_adlo_packunpack/'
dir_out = dirrun_pup + 'xx_hires/'

# read in high-res grid and low-res grid
xc_lr = rdmds(griddir_lr + 'XC')
yc_lr = rdmds(griddir_lr + 'YC')

xc_hr = rdmds(griddir_hr + 'XC')
yc_hr = rdmds(griddir_hr + 'YC')

## NEW BATHY
#hfacc_hr = rdmds(griddir_hr + 'hFacC')
hfacc_hr = rdmds('/scratch/shoshi/labsea_MG_12/grid_hires_cleanbathy/hFacC')

# use linear interpolation and inpaint nans for edges

from skimage.restoration import inpaint

def linear_interp(xx_lr, etan=False):

    points = np.array([xc_lr.ravel(), yc_lr.ravel()]).T

    if etan:
        z = 1
    else:
        z = xx_lr.shape[0]

    xx_hr = np.zeros((z, nyh, nxh))
    for i in range(z):

        if etan:
            values = xx_lr.ravel()
        else:
            values = xx_lr[i].ravel()
        values[values == 0] = np.nan

        if np.count_nonzero(~np.isnan(values)) == 0:
            xx_hr[i] = np.nan
        else:

            if np.count_nonzero(~np.isnan(values)) < 4:
                method = 'nearest'
            else:
                method = 'linear'

            # Create a mask to filter out NaN values in the low-resolution data
            mask = ~np.isnan(values)

            # Interpolate onto the high-resolution grid using linear interpolation
            tmp = griddata(points[mask], values[mask], (xc_hr, yc_hr), method=method)

            # Identify where the high-resolution grid still has NaN values after linear interpolation
            nan_mask = np.isnan(tmp)

            if np.any(nan_mask):
                tmp = inpaint.inpaint_biharmonic(tmp, nan_mask, multichannel=False)

            # set floor and ceiling to that of low-res  
            if etan:
                lmin = np.min(xx_lr)
                lmax = np.max(xx_lr)
            else:
                lmin = np.min(xx_lr[i,:,:])
                lmax = np.max(xx_lr[i,:,:])
            tmp = np.clip(tmp, lmin, lmax)

            if xx_lr.shape[0] == nz:
                xx_hr[i] = tmp * hfacc_hr[i]
            else:
                xx_hr[i] = tmp * hfacc_hr[0]
    print('done')

    return xx_hr


def write_float64(fout,fld):
    with open(fout, 'wb') as f:
        np.array(fld, dtype=">f8").tofile(f)


# read in all xx files

import glob
import os

files = glob.glob(os.path.join(dirrun_pup, 'xx_*.000000' + iter + '*'))
xx_files = []
for f in files:
    if len(f.split('.')) == 3:
        if f.split('.')[-1] == 'meta':
            xx_files.append(f[:-5])



for xx_name in xx_files:
    
    xx = rdmds(xx_name)
    if 'etan' in xx_name:
        xx_hr = linear_interp(xx, etan=True)[0]
    else:
        xx_hr = linear_interp(xx)
    xx_hr[np.isnan(xx_hr)] = 0
    
    write_float64(dir_out + xx_name.split('/')[-1] + '.data', xx_hr)



