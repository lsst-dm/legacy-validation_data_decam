setup obs_decam -t w_2018_37
# To be run from within validation_data_decam repo
# were we were at DM-14868
setup -k -r .  

mkdir data
echo lsst.obs.decam.DecamMapper > data/_mapper
ingestImagesDecam.py data ${VALIDATION_DATA_DECAM_DIR}/instcal/*.fz  --mode copy --filetype instcal
# Link in the reference catalogs
ln -s ${VALIDATION_DATA_DECAM_DIR}/ref_cats data/ref_cats

export OMP_NUM_THREADS=1  # Suppress OMP parallelism.  We parallelize by CCD.
processCcd.py data --output data \
    @${VALIDATION_DATA_DECAM_DIR}/Decam.list \
    --configfile ${VALIDATION_DATA_DECAM_DIR}/decamConfig.py \
    -j 4 \
   >& processCcd.log
