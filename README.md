validation_data_decam
=====================

Data from DECam to validate the performance of the LSST DM stack.
Data from COSMOS field (Dey, NOAO PropID: 2013A-0351)

Data are used by `validate_drp`.

Currently just includes calibrated image data with weightmaps and mask:
Files
------
path      | description
:---------|:-----------------------------
`instcal` | Photometrically and astrometrically calibrated data
          |   as processed by the NOAO DECam Community Pipeline
`dqmask`  | Mask for `instcal` images
`wtmap`   | Weightmap for `instcal` images

No astrometric catalogs are included.

Future directions include
1. Add astrometric and photometric catalogs to re-check the `instcal` calibration
2. Add raw versions of images + calibration data (bias, flat, dark) to enable full ISR processing of data.

Git LFS
-------

To clone and use this repository, you'll need Git Large File Storage (LFS).

Our [Developer Guide](http://developer.lsst.io/en/latest/tools/git_lfs.html) explains how to setup Git LFS for LSST development.
