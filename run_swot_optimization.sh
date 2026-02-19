#!/bin/bash -x
#SBATCH -J labseaMGswot
#SBATCH -o labseaMGswot.%j.out
#SBATCH -e labseaMGswot.%j.err
#SBATCH -t 48:00:00
#SBATCH -p skx
#SBATCH -N 6 
#SBATCH -n 180
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
nprocs_hr=180
nprocs_lr=16

iter=0
itermax=10
costfactor=0.95

jobfile=run_swot_optimization.bash

#--- set dir ------------
#rootdir=/home/shoshi/MITgcm_c69j/lab_sea12/
#scratchdir=/scratch/shoshi/labsea_MG_12/assim_argo_MG
rootdir=/work2/08382/shoshi/stampede3/MITgcm_c69j/lab_sea12/
scratchdir=/scratch/08382/shoshi/labsea_runs/assim_swot_MG/

builddir_hi=${rootdir}/build_adhi_2lev_seaice_update_mpi
builddir_lo=${rootdir}/build_adlo_2lev_seaice_update_mpi
optimdir=${scratchdir}/OPTIM

# --- optim ---

while [ ! ${iter} -gt $itermax ]; do

  ext2=$(printf "%04d" $iter)

# --- 1.high-res forward run ---  
  
  workdir_hi=${scratchdir}/run_adhi_it${ext2}

  if [ ! -d $workdir_hi ]; then
    mkdir -p $workdir_hi;
  fi
  
  cd $workdir_hi
  
  # cp binaries into workdir_hi
  # change data.ctrl
  # cp xx_[ctrl] adjustments if iter > 0
  ${rootdir}/link_hires.sh $iter $ext2 $scratchdir $builddir_hi 

  #---  run  --------
  \rm -f mitgcmuv*
  cp -f ${builddir_hi}/mitgcmuv_ad ./
  cp -f ${builddir_hi}/Makefile ./
  
  set -x
  date > run.MITGCM.timing
#  mpiexec -np ${nprocs_hr} ./mitgcmuv_ad > stdout
  ibrun -n ${nprocs_hr} ./mitgcmuv_ad > stdout
#  ibrun -npernode 16 ./mitgcmuv_ad > stdout
  date >> run.MITGCM.timing
  cd ..

# --- 2.low-res adjoint run ---
  
  workdir_lo=${scratchdir}/run_adlo_it${ext2}

  if [ ! -d $workdir_lo ]; then
    mkdir -p $workdir_lo;
  fi

  source $(conda info --base)/etc/profile.d/conda.sh
  conda activate py38
  # create low-res xx_[ctrl]
  python3 ${rootdir}/mappings/make_cost_cp_v2.py "$ext2" "" "$scratchdir" 
  ## profiles retiling 
  python3 ${rootdir}/mappings/make_obsfit_lr_tiles.py "$ext2" "" 
  conda deactivate

  cd $workdir_lo

  # cp binaries into workdir_lo
  # cp ONLINE low-res cost, misfit, barfiles, and xx_[ctrl]
  # create data.optim
  ${rootdir}/link_lores.sh $iter $workdir_hi $scratchdir $builddir_lo
  
  #---  run  --------
  \rm -f mitgcmuv*
  cp -f ${builddir_lo}/mitgcmuv_ad ./
  cp -f ${builddir_lo}/Makefile ./
  
  set -x
  date > run.MITGCM.timing
#  mpiexec -np ${nprocs_lr} ./mitgcmuv_ad > stdout
  ibrun -n ${nprocs_lr} ./mitgcmuv_ad > stdout
  date >> run.MITGCM.timing

#  sed -i 's/61/0/g' divided.ctrl
  sed -i 's/376/0/g' divided.ctrl
  date > run.MITGCM.timing
#  mpiexec -np ${nprocs_lr} ./mitgcmuv_ad > stdout
  ibrun -n ${nprocs_lr} ./mitgcmuv_ad > stdout
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

# --- 4. pack unpack -----------

  workdir_pup=${scratchdir}/run_adlo_packunpack

  if [ ! -d $workdir_pup ]; then
    mkdir -p $workdir_pup;
  fi

  cd $workdir_pup

  # cp ecco*ctrl*iter from optim 
  # change data.ctrl to unpack=TRUE
  # change data to run for 1 timestep
  # use lo build until figure out why we're not getting xx_[].effective
  # but it's only running the first chunk of the divided adjoint before automatically stopping 
  ${rootdir}/link_packunpack.sh $ext2 $workdir_pup $optimdir $builddir_lo $workdir_hi 

  #---  run  --------
  \rm -f mitgcmuv*
  cp -f ${builddir_lo}/mitgcmuv_ad ./
  cp -f ${builddir_lo}/Makefile ./
  
  set -x
  date > run.MITGCM.timing
#  mpiexec -np ${nprocs_lr} ./mitgcmuv_ad > stdout #2>&1 &
  ibrun -n ${nprocs_lr} ./mitgcmuv_ad > stdout #2>&1 &

#  # Get the PID of the executable
#  EXEC_PID=$!
#  
#  # Monitor the log file in a loop using grep
#  while true; do
#      if grep -q "time_tsnumber" STDOUT.0000; then
#          echo "Pattern found, killing executable."
#          kill $EXEC_PID
#          break
#      fi
#      sleep 1  # Wait for 1 second before checking again
#  done
#  
#  # Optionally, wait for the executable to finish if not killed
#  wait $EXEC_PID

  date >> run.MITGCM.timing
  
  cd ..

# --- 5. interpolate adjustments to high-res -----------
echo $(printf "%04d" $((iter+1)))  # "000$((iter+1))"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate py38
python3 ${rootdir}/mappings/interp_xx_lores_to_hires_itX_v2.py $(printf "%04d" $((iter+1))) $workdir_pup  #"000$((iter+1))" $interp_type
conda deactivate

#  let iter=6
  let iter=iter+1
done

echo "DONE"
