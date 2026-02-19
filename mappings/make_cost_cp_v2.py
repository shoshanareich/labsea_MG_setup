import numpy as np
import sys
sys.path.append('/work2/08382/shoshi/stampede3/MITgcm_c69j/MITgcm/utils/python/MITgcmutils')
#sys.path.append('/home/shoshi/MITgcm_c68r/MITgcm/utils/python/MITgcmutils')
from MITgcmutils import rdmds
from file_utils import *

nx=120
ny=96
nz=21

factor = 4 # lowres * factor = hires

nxh=nx*factor
nyh=ny*factor

#dirroot='/scratch/shoshi/labsea_MG_12/'
dirroot='/scratch/08382/shoshi/labsea_runs/'

#iter = '0000'
iter = sys.argv[1] # sys.argv[0] is name of python file
print(iter)
ext = sys.argv[2]
rundirs = sys.argv[3]

dirrun_lr = rundirs + '/run_adlo_it' + iter + ext + '/' #'_5day_ASTEwts/'
dirrun_hr = rundirs + '/run_adhi_it' + iter + ext + '/' #'_5day_ASTEwts/'


## NEW BATHY
# hfacc_hr = rdmds(dirrun_hr + 'hFacC')
maskc_lr = rdmds(dirroot + 'grid_lores_cleanbathy_v2/maskCtrlC')

def bin_avg(fld, etan=False):
    if etan:
        fld = fld.reshape(ny, factor, nx, factor)
        fld_lr = np.nanmean(fld, axis=(1,3))
    else:
        fld = fld.reshape(fld.shape[0], ny, factor, nx, factor)
        fld_lr = np.nanmean(fld, axis=(2,4))
    return fld_lr

def write_float64(fout,fld):
    with open(fout, 'wb') as f:
        np.array(fld, dtype=">f8").tofile(f)


# read in high-res xx, adxx, m_, adm_

xx_files = get_files(dirrun_hr, 'xx')
adxx_files = get_files(dirrun_hr, 'adxx')
m_files = get_files(dirrun_hr, 'm')
adm_files = get_files(dirrun_hr, 'adm')

# bin-avg to low-res
for f in xx_files.union(xx_files, adxx_files, m_files, adm_files):
    print(f)
    xx = rdmds(dirrun_hr + f)
    xx[xx == 0] = np.nan
    
    etan = 'etan' in f
    if xx.shape[0] == nz:
        xx_lr = bin_avg(xx, etan=etan) * maskc_lr
    else:
        xx_lr = bin_avg(xx, etan=etan) * maskc_lr[0]

    xx_lr[np.isnan(xx_lr)] = 0
    write_float64(dirrun_lr + f + '.data', xx_lr)




#
#####

#misfit_sst = rdmds(dirrun_hr + 'misfit_sst')
#
#misfit_sst[misfit_sst == 0] = np.nan
#
#misfit_sst_lr = bin_avg(misfit_sst) * maskc_lr[0]
#
#misfit_sst_lr[np.isnan(misfit_sst_lr)] = 0
#
#
#write_float64(dirrun_lr + 'misfit_sst.data', misfit_sst_lr)
#
#
#
### confirm high-res cost
#sigma = 0.010402
#cost_hr = np.nansum((misfit_sst / sigma)**2)
#print(cost_hr)
#
#cost_lr = np.nansum((misfit_sst_lr * maskc_lr[0] / sigma) **2)
#print(f'total cost = {cost_lr}')
#
#with open(dirrun_hr + 'costfinal_lo', 'w') as f:
#    f.write(str(cost_lr))
