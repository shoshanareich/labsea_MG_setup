import numpy as np
import sys
sys.path.append('/home/shoshi/MITgcm_c68r/MITgcm/utils/python/MITgcmutils')
from MITgcmutils import rdmds

nx=120
ny=96
nz=21

factor = 4 # lowres * factor = hires

nxh=nx*factor
nyh=ny*factor

dirroot='/scratch/shoshi/labsea_MG_12/'

#iter = '0000'
iter = sys.argv[1] # sys.argv[0] is name of python file
print(iter)

dirrun_lr = dirroot + 'assim_swot_argo/run_adlo_it' + iter + '/'
dirrun_hr = dirroot + 'assim_swot_argo/run_adhi_it' + iter + '/'


## NEW BATHY
# hfacc_hr = rdmds(dirrun_hr + 'hFacC')
maskc_lr = rdmds(dirroot + 'grid_lores/maskCtrlC')

def bin_avg(fld):
    fld = fld.reshape(fld.shape[0], ny, factor, nx, factor)
    fld_lr = np.nanmean(fld, axis=(2,4))
    return fld_lr

def write_float64(fout,fld):
    with open(fout, 'wb') as f:
        np.array(fld, dtype=">f8").tofile(f)


# read in high-res 
xx_theta = rdmds(dirrun_hr + 'xx_theta.000000' + iter)
xx_atemp = rdmds(dirrun_hr + 'xx_atemp.000000' + iter)

adxx_theta = rdmds(dirrun_hr + 'adxx_theta.000000' + iter)
adxx_atemp = rdmds(dirrun_hr + 'adxx_atemp.000000' + iter)

m_sst_day = rdmds(dirrun_hr + 'm_sst_day.000000' + iter)
adm_sst_day = rdmds(dirrun_hr + 'adm_sst_day.000000' + iter)

misfit_sst = rdmds(dirrun_hr + 'misfit_sst')



# bin-avg to low-res
xx_theta[xx_theta == 0] = np.nan
xx_atemp[xx_atemp == 0] = np.nan

adxx_theta[adxx_theta == 0] = np.nan
adxx_atemp[adxx_atemp == 0] = np.nan

m_sst_day[m_sst_day == 0] = np.nan
adm_sst_day[adm_sst_day == 0] = np.nan

misfit_sst[misfit_sst == 0] = np.nan

xx_theta_lr = bin_avg(xx_theta) * maskc_lr 
xx_atemp_lr = bin_avg(xx_atemp) * maskc_lr[0]

adxx_theta_lr = bin_avg(adxx_theta) * maskc_lr
adxx_atemp_lr = bin_avg(adxx_atemp) * maskc_lr[0]

m_sst_day_lr = bin_avg(m_sst_day) * maskc_lr[0]
adm_sst_day_lr = bin_avg(adm_sst_day) * maskc_lr[0]

misfit_sst_lr = bin_avg(misfit_sst) * maskc_lr[0]


xx_theta_lr[np.isnan(xx_theta_lr)] = 0
xx_atemp_lr[np.isnan(xx_atemp_lr)] = 0

adxx_theta_lr[np.isnan(adxx_theta_lr)] = 0
adxx_atemp_lr[np.isnan(adxx_atemp_lr)] = 0

m_sst_day_lr[np.isnan(m_sst_day_lr)] = 0
adm_sst_day_lr[np.isnan(adm_sst_day_lr)] = 0

misfit_sst_lr[np.isnan(misfit_sst_lr)] = 0


write_float64(dirrun_lr + 'xx_theta.000000' + iter + '.data', xx_theta_lr)
write_float64(dirrun_lr + 'xx_atemp.000000' + iter + '.data', xx_atemp_lr)
write_float64(dirrun_lr + 'adxx_theta.000000' + iter + '.data', adxx_theta_lr)
write_float64(dirrun_lr + 'adxx_atemp.000000' + iter + '.data', adxx_atemp_lr)
write_float64(dirrun_lr + 'm_sst_day.000000' + iter + '.data', m_sst_day_lr)
write_float64(dirrun_lr + 'adm_sst_day.000000' + iter + '.data', adm_sst_day_lr)
write_float64(dirrun_lr + 'misfit_sst.data', misfit_sst_lr)



# confirm high-res cost
sigma = 0.010402
cost_hr = np.nansum((misfit_sst / sigma)**2)
print(cost_hr)

cost_lr = np.nansum((misfit_sst_lr * maskc_lr[0] / sigma) **2)
print(f'total cost = {cost_lr}')

with open(dirrun_hr + 'costfinal_lo', 'w') as f:
    f.write(str(cost_lr))
