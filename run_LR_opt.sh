#!/bin/bash -x
#SBATCH -J labseaLRargo 
#SBATCH -o labseaLRargo.%j.out
#SBATCH -e labseaLRargo.%j.err
#SBATCH -t 6:00:00
#SBATCH -p skx
#SBATCH -N 1 
#SBATCH -n 16
#SBATCH -A OCE23001
#SBATCH --mail-user=sreich@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

## SVERDRUP
##SBATCH -N 8
##SBATCH -n 180

#--- 0.load modules ------
#module purge
##module load intel openmpi netcdf-fortran
#module load intel/2023.1.0 openmpi4/4.1.5 phdf5/1.14.1 netcdf-fortran/4.6.0 netcdf/4.9.0 prun
#echo $LD_LIBRARY_PATH

module purge; module load intel/25.1 impi/21.15 netcdf/4.9.2 hdf5/1.14.6

#ulimit -s hard
#ulimit -u hard
ulimit -s unlimited
ulimit -v unlimited
#export IBRUN_TASKS_PER_NODE=16
export I_MPI_DEBUG=4

#export UCX_MEMTYPE_CACHE=n
#export UCX_TLS=rc,self,sm

#---- set variables ------
# note: for nprocs, take ntiles - length(blanklist)
nprocs_lr=16

iter=0
itermax=10
costfactor=0.95

jobfile=run_optimization.bash

#--- set dir ------------
#rootdir=/home/shoshi/MITgcm_c69j/lab_sea12/
#scratchdir=/scratch/shoshi/labsea_MG_12/assim_argo_MG
rootdir=/work2/08382/shoshi/stampede3/MITgcm_c69j/lab_sea12/
scratchdir=/scratch/08382/shoshi/labsea_runs/assim_argo_LR/

builddir_lo=${rootdir}/build_adlo_2lev_seaice_update_mpi
optimdir=${scratchdir}/OPTIM

# --- optim ---

while [ ! ${iter} -gt $itermax ]; do

  ext2=$(printf "%04d" $iter)

# --- low-res forward and adjoint run ---
  
  workdir_lo=${scratchdir}/run_lo_it${ext2}

  if [ ! -d $workdir_lo ]; then
    mkdir -p $workdir_lo;
  fi

  cd $workdir_lo

  mkdir -p ./diags/     

  #--- 6. NAMELISTS ---------
  #ln -s ${rootdir}/input_cal/* .
  cp ${rootdir}/input_adlo/* .
  ln -s ${rootdir}/input_binaries_lores/bathy_cleaned_96x120.bin .
  ln -s ${rootdir}/input_binaries_lores/LevCli_temp_120x96_linearv3_smooth.labsea1979 .
  ln -s ${rootdir}/input_binaries_lores/LevCli_salt_120x96_linearv3_smooth.labsea1979 .
  ln -s ${rootdir}/input_binaries_lores/viscd_bottom_topo*.bin .
  ln -s ${rootdir}/input_binaries_lores/viscz_bottom_topo*.bin .
  ln -s ${rootdir}/input_binaries_lores/diffkr_r4.bin .
  ln -s ${rootdir}/input_binaries_exf/* .
  ln -s ${rootdir}/input_binaries_hires/ones_64b.bin .
  ln -s ${rootdir}/input_binaries_hires/ARGO_WO_2024_PFL_D_labsea_splitcost.nc .
  ln -s ${rootdir}/input_weights_lores/*_jra3q_weights_Jan2024_64b_SMOOTHED_removeboundary.bin .
  ln -s ${rootdir}/input_binaries_lores/rads_j3_labsea_96x120_v2_2024 .
  ln -s ${rootdir}/input_binaries_lores/slaerr_03m.bin .
  ln -s ${rootdir}/input_weights_lores/*fromASTE_*.bin .
  cp ${rootdir}/input_binaries_lores/pickup_seaice.0000025920.meta pickup_seaice.0000000001.meta
  cp ${rootdir}/input_binaries_lores/pickup_seaice.0000025920.data pickup_seaice.0000000001.data
  cp ${rootdir}/input_binaries_lores/pickup.0000025920.meta pickup.0000000001.meta
  cp ${rootdir}/input_binaries_lores/pickup.0000025920.data pickup.0000000001.data
  
  mkdir jra3q
  #ln -s /scratch/shared/jra3q/jra3q_*_2024 ./jra3q/
  ln -s /scratch/08382/shoshi/jra3q/jra3q_*_2024 ./jra3q/
  rm data.diagnostics
  cp ${rootdir}/input_adhi/data.diagnostics .
  
  ##--- 5. linking xx_ fields ------
  if [ ${iter} -lt 1 ]; then
    sed -i -e 's/'"doinitxx = .FALSE."'/'"doinitxx = .TRUE."'/g' data.ctrl
    sed -i -e 's/'"doInitXX = .FALSE."'/'"doInitXX = .TRUE."'/g' data.ctrl
    sed -i -e 's/'"doMainPack = .FALSE."'/'"doMainPack = .TRUE."'/g' data.ctrl
    sed -i -e 's/'"doMainUnpack = .TRUE."'/'"doMainUnpack = .FALSE."'/g' data.ctrl
  else
    sed -i -e 's/'"doinitxx = .TRUE."'/'"doinitxx = .FALSE."'/g' data.ctrl
    sed -i -e 's/'"doInitXX = .TRUE."'/'"doInitXX = .FALSE."'/g' data.ctrl
    sed -i -e 's/'"doMainPack = .FALSE."'/'"doMainPack = .TRUE."'/g' data.ctrl
    sed -i -e 's/'"doMainUnpack = .FALSE."'/'"doMainUnpack = .TRUE."'/g' data.ctrl
    cp ${optimdir}/ecco_ctrl_MIT_CE_000.opt${ext2} .
  fi
  #--- 10. (re)set optimcycle --------------------
  
  \rm data.optim
  cat > data.optim <<EOF
   &OPTIM
   optimcycle=${iter},
   /
EOF

  #---  run forward  --------
  cp -f ${builddir_lo}/mitgcmuv_ad ./
  cp -f ${builddir_lo}/Makefile ./
  
  set -x
  date > run.MITGCM.timing
#  mpiexec -np ${nprocs_lr} ./mitgcmuv_ad > stdout
  ibrun -n ${nprocs_lr} ./mitgcmuv_ad > stdout # run forward
  ibrun -n ${nprocs_lr} ./mitgcmuv_ad > stdout # run first chunk of DivA

  sed -i 's/376/0/g' divided.ctrl 
  ibrun -n ${nprocs_lr} ./mitgcmuv_ad > stdout # run rest of DivA
  date >> run.MITGCM.timing
  cd ..


# --- 3. OPTIM ----------------
  cd ${optimdir}
  #bash reset.bash
  cp ${workdir_lo}/ecco_cost_MIT_CE_000.opt${ext2} ${optimdir} 
  cp ${workdir_lo}/ecco_ctrl_MIT_CE_000.opt${ext2} ${optimdir} 
  cp ${workdir_lo}/data.ctrl ${optimdir} 
  cp -f ${workdir_lo}/costfinal ${optimdir}
#  cost=$(grep fc costfunction${ext2}  | sed 's/D/E/g' | awk '{printf "%14.12e", $3}')
#  costf=$(grep fc costfunction${ext2} | sed 's/D/E/g' | awk '{printf "%0.14f", $3}')
  echo "iter = $iter"
#  echo "cost = $cost"
  
#  costupdate=$(echo $costf*$costfactor | bc)
#  costnew=`echo $costupdate|awk '{printf "%14.12e\n", $costupdate}'`
#  echo "costnew = $costnew"

  mv data.optim data.optim_bk
  cat > data.optim <<EOF
    &OPTIM
    optimcycle=${iter},
    numiter=100,
    nfunc=9,
    dfminFrac=0.05,
#    fmin=${costnew},
    iprint=10,
    nupdate=4,
    /
    &M1QN3
    coldstart = .TRUE.,
    /
EOF

  \rm OP*
  ./optim.x > opt_it${ext2}.txt

  dir_iter=${workdir_lo}/optim/+it${ext2}
  mkdir -p $dir_iter
  cp data.optim $dir_iter
  cp data.optim data.optim_${ext2}
  cp ecco_ctrl_MIT_CE_000.opt${ext2} $dir_iter

  cd ..

  let iter=iter+1
done

echo "DONE"
