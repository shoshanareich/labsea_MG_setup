rootdir=/home/shoshi/MITgcm_c69j/lab_sea12/
datadir=/home/shoshi/MITgcm_obsfit/lab_sea12/
scratchdir=/scratch/shoshi/labsea_MG_12/assim_argo_MG/
iter=$1
iter=$((iter+1))
optimext=$(printf "%04d" $iter)
workdir=$2
optimdir=$3
builddir=$4
dirhires=$5

mkdir -p ./diags/

#--- 6. NAMELISTS ---------
#ln -s ${datadir}/input_cal/* .
cp ${rootdir}/input_adlo/* .
ln -s ${datadir}/input_binaries_lores/bathy_cleaned_96x120.bin .
ln -s ${datadir}/input_binaries_lores/LevCli_temp_120x96_linearv3_smooth.labsea1979 .
ln -s ${datadir}/input_binaries_lores/LevCli_salt_120x96_linearv3_smooth.labsea1979 .
ln -s ${datadir}/input_binaries_lores/viscd_bottom_topo*.bin .
ln -s ${datadir}/input_binaries_lores/viscz_bottom_topo*.bin .
ln -s ${datadir}/input_binaries_lores/diffkr_r4.bin .
ln -s ${datadir}/input_binaries_exf/* .
ln -s ${datadir}/input_binaries_hires/ones_64b.bin .
ln -s ${datadir}/input_binaries_hires/ARGO_WO_2024_PFL_D_labsea_splitcost.nc .
ln -s ${datadir}/input_weights_lores/*_jra3q_weights_Jan2024_64b_SMOOTHED_removeboundary.bin .
ln -s ${datadir}/input_binaries_lores/rads_j3_labsea_96x120_v2_2024 .
ln -s ${datadir}/input_binaries_lores/slaerr_03m.bin .
ln -s ${datadir}/input_weights_lores/*fromASTE_*.bin .
cp /scratch/shoshi/labsea_MG_12/run_lo_1yr_jra3q_B/pickup_seaice.0000025920.meta pickup_seaice.0000000001.meta
cp /scratch/shoshi/labsea_MG_12/run_lo_1yr_jra3q_B/pickup_seaice.0000025920.data pickup_seaice.0000000001.data
cp /scratch/shoshi/labsea_MG_12/run_lo_1yr_jra3q_B/pickup.0000025920.meta pickup.0000000001.meta
cp /scratch/shoshi/labsea_MG_12/run_lo_1yr_jra3q_B/pickup.0000025920.data pickup.0000000001.data

mkdir jra3q
ln -s /scratch/shared/jra3q/jra3q_*_2024 ./jra3q/
rm data.diagnostics
cp ${datadir}/input_adhi/data.diagnostics .

cp ${optimdir}/ecco_ctrl_MIT_CE_000.opt$optimext .
rm costfinal

#tapes:
rm -rf ./tapes
cp -r $dirhires/tapes ./

##--- 5. linking xx_ fields ------
    sed -i -e 's/'"doinitxx = .TRUE"'/'"doinitxx = .FALSE"'/g' data.ctrl
    sed -i -e 's/'"doInitxx = .TRUE"'/'"doInitxx = .FALSE"'/g' data.ctrl
    sed -i -e 's/'"doInitXX = .TRUE."'/'"doInitXX = .FALSE."'/g' data.ctrl
    sed -i -e 's/'"doMainUnpack = .FALSE."'/'"doMainUnpack = .TRUE."'/g' data.ctrl
    sed -i -e 's/'"doMainUnpack = .false."'/'"doMainUnpack = .true."'/g' data.ctrl
    sed -i -e 's/'"doMainPack = .TRUE."'/'"doMainPack = .FALSE."'/g' data.ctrl
    sed -i -e 's/'"doMainPack = .true."'/'"doMainPack = .false."'/g' data.ctrl

#--- 10. (re)set optimcycle --------------------

\rm data.optim
cat > data.optim <<EOF
 &OPTIM
 optimcycle=${iter},
 /
EOF

#--- 7. executable --------
cp -p ${builddir}/mitgcmuv_ad .

mkdir xx_hires/
