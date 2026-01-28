iter=$1
optimext=$2

#--- 2.set dir ------------
rootdir=/home/shoshi/MITgcm_c69j/lab_sea12/
datadir=/home/shoshi/MITgcm_obsfit/lab_sea12/
scratchdir=/scratch/shoshi/labsea_MG_12/assim_argo_MG/
dirrun_pup=${scratchdir}/run_adlo_packunpack/

hr_dir_iter0=${scratchdir}/run_adhi_it0000

builddir=$3

mkdir diags

#--- 6. NAMELISTS ---------
#ln -s ${datadir}/input_cal/* .
cp ${rootdir}/input_adhi/* .
ln -s ${datadir}/input_binaries_hires/bathy_cleaned.bin .
ln -s ${datadir}/input_binaries_hires/LevCli_temp_480x384_linearv5.labsea1979 .
ln -s ${datadir}/input_binaries_hires/LevCli_salt_480x384_linearv5.labsea1979 .
ln -s ${datadir}/input_binaries_exf/* .
ln -s ${datadir}/input_binaries_hires/ones_64b.bin .
ln -s ${datadir}/input_binaries_hires/ARGO_WO_2024_PFL_D_labsea_splitcost.nc .
ln -s ${datadir}/input_weights_hires/*_jra3q_weights_Jan2024_64b_SMOOTHED_removeboundary.bin .
ln -s ${datadir}/input_weights_hires/*fromASTE_*.bin .
ln -s ${datadir}/input_binaries_hires/diffkr_r4_HR.bin .
cp /scratch/shoshi/labsea_MG_12//run_hi_1yr_jra3q_B/pickup_seaice.0000103680.meta pickup_seaice.0000000001.meta
cp /scratch/shoshi/labsea_MG_12//run_hi_1yr_jra3q_B/pickup_seaice.0000103680.data pickup_seaice.0000000001.data
cp /scratch/shoshi/labsea_MG_12//run_hi_1yr_jra3q_B/pickup.0000103680.meta pickup.0000000001.meta
cp /scratch/shoshi/labsea_MG_12//run_hi_1yr_jra3q_B/pickup.0000103680.data pickup.0000000001.data
cp ${datadir}/input_binaries_hires/smooth2Dscales001_4x.bin ./smooth2Dscales000.data
cp ${datadir}/input_binaries_hires/smooth2Dscales001_4x.bin ./smooth2Dscales001.data
cp ${datadir}/input_binaries_hires/smooth3DscalesH001_4x.bin ./smooth3DscalesH001.data
cp ${datadir}/input_binaries_hires/smooth3DscalesZ001_4x.bin ./smooth3DscalesZ001
cp ${datadir}/input_binaries_hires/smooth*operator* .
cp ${datadir}/input_binaries_hires/smooth2Doperator001.meta ./smooth2Doperator000.meta
cp ${datadir}/input_binaries_hires/smooth2Doperator001.data ./smooth2Doperator000.data
cp ${datadir}/input_binaries_hires/smooth2Dnorm001.meta ./smooth2Dnorm000.meta
cp ${datadir}/input_binaries_hires/smooth2Dnorm001.data ./smooth2Dnorm000.data
cp ${datadir}/input_binaries_hires/smooth*norm* .
mkdir jra3q
ln -s /scratch/shared/jra3q/jra3q_*_2024 ./jra3q/

#-- swap out data.ctrl and copy high-res adjustments
if [ ${iter} -lt 1 ]; then
  sed -i -e 's/'"doinitxx = .FALSE."'/'"doinitxx = .TRUE."'/g' data.ctrl
  sed -i -e 's/'"doInitXX = .FALSE."'/'"doInitXX = .TRUE."'/g' data.ctrl
  sed -i -e 's/'"doMainPack = .FALSE."'/'"doMainPack = .TRUE."'/g' data.ctrl
  sed -i -e 's/'"doMainUnpack = .TRUE."'/'"doMainUnpack = .FALSE."'/g' data.ctrl
else
#  mkdir adxxfiles/
  cp ${dirrun_pup}/xx_hires/*.data .
  cp -f ${datadir}/input_binaries_hires/xx_theta.0000000000.meta ./xx_theta.000000${optimext}.meta
  cp -f ${datadir}/input_binaries_hires/xx_atemp.0000000000.meta ./xx_atemp.000000${optimext}.meta
  cp $hr_dir_iter0/xx_theta.0000000000.meta ./xx_theta.000000${optimext}.meta
  cp $hr_dir_iter0/xx_salt.0000000000.meta ./xx_salt.000000${optimext}.meta
  cp $hr_dir_iter0/xx_atemp.0000000000.meta ./xx_atemp.000000${optimext}.meta
  cp $hr_dir_iter0/xx_precip.0000000000.meta ./xx_precip.000000${optimext}.meta
  cp $hr_dir_iter0/xx_swdown.0000000000.meta ./xx_swdown.000000${optimext}.meta
  cp $hr_dir_iter0/xx_lwdown.0000000000.meta ./xx_lwdown.000000${optimext}.meta
  cp $hr_dir_iter0/xx_uwind.0000000000.meta ./xx_uwind.000000${optimext}.meta
  cp $hr_dir_iter0/xx_vwind.0000000000.meta ./xx_vwind.000000${optimext}.meta
  cp $hr_dir_iter0/xx_aqh.0000000000.meta ./xx_aqh.000000${optimext}.meta
  cp $hr_dir_iter0/xx_diffkr.0000000000.meta ./xx_diffkr.000000${optimext}.meta
  sed -i -e 's/'"doinitxx = .TRUE."'/'"doinitxx = .FALSE."'/g' data.ctrl
  sed -i -e 's/'"doInitXX = .TRUE."'/'"doInitXX = .FALSE."'/g' data.ctrl
  sed -i -e 's/'"doMainPack = .FALSE."'/'"doMainPack = .TRUE."'/g' data.ctrl
  sed -i -e 's/'"doMainUnpack = .FALSE."'/'"doMainUnpack = .FALSE."'/g' data.ctrl
  sed -i -e 's/'"doMainUnpack = .false."'/'"doMainUnpack = .false."'/g' data.ctrl

\rm data.optim
cat > data.optim <<EOF
 &OPTIM
 optimcycle=${iter},
 /
EOF

fi

#--- 7. executable --------
cp -p $builddir/mitgcmuv_ad .


