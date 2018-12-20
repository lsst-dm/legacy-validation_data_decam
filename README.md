validation_data_decam
=====================

Data from DECam to validate the performance of the LSST DM stack.
Data from COSMOS field (Dey, NOAO PropID: 2013A-0351)

Data are targeted at use by LSST DM Stack developers, continuous integration tests, and `validate_drp`.  General users who wish to simply see examples of the steps and products of the LSST DM Stack processing may likely also find this useful.

The use model for reprocessing these data for comparisons and continuous integration, is to use these data in as an input repository and then write to a separate output repository (by, e.g.,  specifying `--output` on the command line call to a Task).  One might then wish to compare the outputs of the two repositories at the image, catalog, or summary statistic or visualization plot level.

These data are provided as a Butler repository in its POSIX filesystem mode*.  This provides straightforward access through both the Butler mechanism for clean programatic use, but also relatively simple access directly from the file system for quick inspection and orientation to the data.

[*] Mode isn't quite the right word here.

The DECam data are currently just "Community Pipeline" calibrated image data with weightmaps and mask:

Files
------
path                  | description
:---------------------|:-----------------------------
`instcal`             | Photometrically and astrometrically calibrated data
                      |   as processed by the NOAO DECam Community Pipeline
`dqmask`              | Mask for `instcal` images
`wtmap`               | Weightmap for `instcal` images
`data`                | Butler repo of ingested raw data and processCcd results
`astrometry_net_data` | SDSS DR9 catalog files in astrometry.net format
                      |   as photometrically recalibrated by Doug Finkbeiner's group.
`ref_cats`            | HTM indexed catalog files from both SDSS and Pan-Starrs for
                      |   astrometric and photometric calibration.
`processCcd.log`      | Output from the processCcd run on the `data` repo.
`eups_setup_used.txt` | EUPS setup configuration for ingestImages and processCcd run
`Decam.list`          | List of dataIds in this repo.  For use in running Tasks.

No astrometric catalogs are included.

Future directions include
1. Add astrometric and photometric catalogs to re-check the `instcal` calibration
2. Add raw versions of images + calibration data (bias, flat, dark) to enable full ISR processing of data.

Example usage
-------------

```
setup validation_data_decam

NEW_OUTPUT_REPO=DECam_data

processCcd.py ${VALIDATION_DATA_DECAM_DIR}/data \
    --output ${NEW_OUTPUT_REPO} \
    --configfile ${VALIDATION_DATA_DECAM_DIR}/decamConfig.py \
    @${VALIDATION_DATA_DECAM_DIR}/Decam.list \
    -j 4 
```

Notes:
 * There will be a `${NEW_OUTPUT_REPO}/_parent` link back to the input repository `${VALIDATION_DATA_DECAM_DIR}/data`.
 * The list of images (`dataIds`) to process is in `@${VALIDATION_DATA_DECAM_DIR}/Decam.list`
 * `-j 4` specifies using 4 cores.  You may wish to change to an appropriate number on your system, but the intent is that `-j 4` should be a reasonable default in 2016.
 * We specifically use `--clobber-config` here because we're running off an already existing repository.

Analyzing the repository
------------------------
One might then choose to use `validate_drp` to analyze the peformance of the results against SRD metrics.

```
setup validate_drp

validateDrp.py DECam_data
```

Recreating the repository
-------------------------
To fully recreate this Butler `repo` from the `raw` data, run the `reprocess.sh` script in this directory.

Notes
 1. We use `--copy` to create a full copy of the raw images in the repo.
 2. No separate `--output` repo was specified.  In this case we intentionally wish to create
the output products in the same repo.

The packages and versions used were recorded into `eups_setup_used.txt`:

```
eups list --setup | awk '{printf "%-30s %s\n", $1, $2}' > eups_setup_used.txt
```

and the repo was made read-only to prevent accidental writes to this repo.

```
chmod -R ugo-w data
```


Git LFS
-------

To clone and use this repository, you'll need Git Large File Storage (LFS).

Our [Developer Guide](http://developer.lsst.io/en/latest/tools/git_lfs.html) explains how to setup Git LFS for LSST development.
