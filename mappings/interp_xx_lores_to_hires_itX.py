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

dirrun_pup = dirroot + 'run_adlo_packunpack_codeupdate/'

# read in unpacked low-res adjustments
xx_theta = rdmds(dirrun_pup + 'xx_theta.000000' + iter).reshape(nz, ny, nx)
xx_atemp = rdmds(dirrun_pup + 'xx_atemp.000000' + iter).reshape(8,ny, nx)


# read in high-res grid and low-res grid
xc_lr = rdmds(griddir_lr + 'XC')
yc_lr = rdmds(griddir_lr + 'YC')

xc_hr = rdmds(griddir_hr + 'XC')
yc_hr = rdmds(griddir_hr + 'YC')

## NEW BATHY
hfacc_hr = rdmds(griddir_hr + 'hFacC')

# use linear interpolation and inpaint nans for edges

from skimage.restoration import inpaint

def linear_interp(xx_lr, dim):

    points = np.array([xc_lr.ravel(), yc_lr.ravel()]).T
    
    xx_hr = np.zeros((xx_lr.shape[0], nyh, nxh))
    for i in range(xx_lr.shape[0]):
        
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
            lmin = np.min(xx_lr[i,:,:])
            lmax = np.max(xx_lr[i,:,:])
            tmp = np.clip(tmp, lmin, lmax)
        
            if dim == 2:
                xx_hr[i] = tmp * hfacc_hr[0]
            elif dim == 3:
                xx_hr[i] = tmp * hfacc_hr[i]
    
    print('done')

    return xx_hr



def write_float64(fout,fld):
    with open(fout, 'wb') as f:
        np.array(fld, dtype=">f8").tofile(f)


xx_theta_hr = linear_interp(xx_theta, 3)
xx_atemp_hr = linear_interp(xx_atemp, 2)

xx_theta_hr[np.isnan(xx_theta_hr)] = 0
xx_atemp_hr[np.isnan(xx_atemp_hr)] = 0

write_float64(dirrun_pup + 'xx_theta_hr.000000' + iter + '.data', xx_theta_hr)
write_float64(dirrun_pup + 'xx_atemp_hr.000000' + iter + '.data', xx_atemp_hr)

