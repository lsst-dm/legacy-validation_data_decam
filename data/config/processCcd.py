import lsst.pipe.tasks.processCcd
assert type(config)==lsst.pipe.tasks.processCcd.ProcessCcdConfig, 'config is of type %s.%s instead of lsst.pipe.tasks.processCcd.ProcessCcdConfig' % (type(config).__module__, type(config).__name__)
import eups.Product
import eups.lock
import eups.stack
import eups.app
import eups.distrib.Repositories
import eups.distrib.server
import eups.utils
import lsst.obs.decam.isr
import eups.db.ChainFile
import eups.distrib.Distrib
import lsst.obs.decam.decamNullIsr
import eups.tags
import eups.cmd
import eups.table
import eups.hooks
import configparser
import eups.distrib.DistribFactory
import eups.db
import eups.distrib.Repository
import eups.db.VersionFile
import optparse
import eups.VersionCompare
import eups.distrib
import eups.distrib.builder
import pipes
import eups.db.Database
import eups.stack.ProductStack
import eups.distrib.eupspkg
import eups.exceptions
import eups.Eups
import eups.Uses
import eups.VersionParser
import eups.distrib.tarball
import eups.distrib.pacman
import eups.stack.ProductFamily
import eups
import lsst.pipe.tasks.setConfigFromEups
import lsst.obs.decam.decamNullIsr
config.isr.retarget(target=lsst.obs.decam.decamNullIsr.DecamNullIsrTask, ConfigClass=lsst.obs.decam.decamNullIsr.DecamNullIsrConfig)

# Persist loaded data as a postISRCCD? The default is false, to avoid duplicating data.
config.isr.doWrite=False

# Dataset type for input data; read by ProcessCcdTask; users will typically leave this alone
config.isr.datasetType='instcal'

# Measure PSF? If False then for all subsequent operations use either existing PSF model when present, or install simple PSF model when not (see installSimplePsf config options)
config.charImage.doMeasurePsf=True

# Persist results?
config.charImage.doWrite=True

# Write icExp and icExpBackground in addition to icSrc? Ignored if doWrite False.
config.charImage.doWriteExposure=False

# Number of iterations of detect sources, measure sources, estimate PSF. If useSimplePsf is True then 2 should be plenty; otherwise more may be wanted.
# 	Valid Range = [1,inf)
config.charImage.psfIterations=2

# type of statistic to use for grid points
# Allowed values:
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 	MEDIAN	median
# 	None	Field is optional
# 
config.charImage.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 
config.charImage.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.charImage.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
# 	Valid Range = [0,inf)
config.charImage.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
# 	Valid Range = [0,inf)
config.charImage.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NONE	No background estimation is to be attempted
# 	None	Field is optional
# 
config.charImage.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.charImage.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.charImage.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.charImage.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.charImage.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.charImage.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.charImage.background.weighting=True

# detected sources with fewer than the specified number of pixels will be ignored
# 	Valid Range = [0,inf)
config.charImage.detection.minPixels=1

# Pixels should be grown as isotropically as possible (slower)
config.charImage.detection.isotropicGrow=False

# Grow detections by nSigmaToGrow * sigma; if 0 then do not grow
config.charImage.detection.nSigmaToGrow=2.4

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.charImage.detection.returnOriginalFootprints=False

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.charImage.detection.thresholdValue=5.0

# Include threshold relative to thresholdValue
# 	Valid Range = [0.0,inf)
config.charImage.detection.includeThresholdMultiplier=10.0

# specifies the desired flavor of Threshold
# Allowed values:
# 	variance	threshold applied to image variance
# 	stdev	threshold applied to image std deviation
# 	value	threshold applied to image value
# 	pixel_stdev	threshold applied to per-pixel std deviation
# 
config.charImage.detection.thresholdType='stdev'

# specifies whether to detect positive, or negative sources, or both
# Allowed values:
# 	positive	detect only positive sources
# 	negative	detect only negative sources
# 	both	detect both positive and negative sources
# 
config.charImage.detection.thresholdPolarity='positive'

# Fiddle factor to add to the background; debugging only
config.charImage.detection.adjustBackground=0.0

# Estimate the background again after final source detection?
config.charImage.detection.reEstimateBackground=True

# type of statistic to use for grid points
# Allowed values:
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 	MEDIAN	median
# 	None	Field is optional
# 
config.charImage.detection.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 
config.charImage.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.charImage.detection.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
# 	Valid Range = [0,inf)
config.charImage.detection.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
# 	Valid Range = [0,inf)
config.charImage.detection.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NONE	No background estimation is to be attempted
# 	None	Field is optional
# 
config.charImage.detection.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.charImage.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.charImage.detection.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.charImage.detection.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.charImage.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.charImage.detection.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.charImage.detection.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 	MEDIAN	median
# 	None	Field is optional
# 
config.charImage.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 
config.charImage.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.charImage.detection.tempLocalBackground.binSize=64

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
# 	Valid Range = [0,inf)
config.charImage.detection.tempLocalBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
# 	Valid Range = [0,inf)
config.charImage.detection.tempLocalBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NONE	No background estimation is to be attempted
# 	None	Field is optional
# 
config.charImage.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.charImage.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.charImage.detection.tempLocalBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.charImage.detection.tempLocalBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.charImage.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.charImage.detection.tempLocalBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.charImage.detection.tempLocalBackground.weighting=True

# Do temporary interpolated background subtraction before footprint detection?
config.charImage.detection.doTempLocalBackground=False

# The maximum number of peaks in a Footprint before trying to replace its peaks using the temporary local background
config.charImage.detection.nPeaksMaxSimple=1

# Run deblender input exposure
config.charImage.doDeblend=False

# What to do when a peak to be deblended is close to the edge of the image
# Allowed values:
# 	clip	Clip the template at the edge AND the mirror of the edge.
# 	ramp	Ramp down flux at the image edge by the PSF
# 	noclip	Ignore the edge when building the symmetric template.
# 	None	Field is optional
# 
config.charImage.deblend.edgeHandling='ramp'

# When the deblender should attribute stray flux to point sources
# Allowed values:
# 	necessary	When there is not an extended object in the footprint
# 	always	Always
# 	never	Never; stray flux will not be attributed to any deblended child if the deblender thinks all peaks look like point sources
# 	None	Field is optional
# 
config.charImage.deblend.strayFluxToPointSources='necessary'

# Assign stray flux (not claimed by any child in the deblender) to deblend children.
config.charImage.deblend.assignStrayFlux=True

# How to split flux among peaks
# Allowed values:
# 	r-to-peak	~ 1/(1+R^2) to the peak
# 	r-to-footprint	~ 1/(1+R^2) to the closest pixel in the footprint.  CAUTION: this can be computationally expensive on large footprints!
# 	nearest-footprint	Assign 100% to the nearest footprint (using L-1 norm aka Manhattan distance)
# 	trim	Shrink the parent footprint to pixels that are not assigned to children
# 	None	Field is optional
# 
config.charImage.deblend.strayFluxRule='trim'

# When splitting stray flux, clip fractions below this value to zero.
config.charImage.deblend.clipStrayFluxFraction=0.001

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.charImage.deblend.psfChisq1=1.5

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.charImage.deblend.psfChisq2=1.5

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.charImage.deblend.psfChisq2b=1.5

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.charImage.deblend.maxNumberOfPeaks=0

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.charImage.deblend.maxFootprintArea=1000000

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.charImage.deblend.maxFootprintSize=0

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.charImage.deblend.minFootprintAxisRatio=0.0

# Mask name for footprints not deblended, or None
config.charImage.deblend.notDeblendedMask='NOT_DEBLENDED'

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
# 	Valid Range = [2,inf)
config.charImage.deblend.tinyFootprintSize=2

# Guarantee that all peaks produce a child source.
config.charImage.deblend.propagateAllPeaks=False

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.charImage.deblend.catchFailures=False

# Mask planes to ignore when performing statistics
config.charImage.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.charImage.deblend.maskLimits={}

# If true, a least-squares fit of the templates will be done to the full image. The templates will be re-weighted based on this fit.
config.charImage.deblend.weightTemplates=False

# Try to remove similar templates?
config.charImage.deblend.removeDegenerateTemplates=False

# If the dot product between two templates is larger than this value, we consider them to be describing the same object (i.e. they are degenerate).  If one of the objects has been labeled as a PSF it will be removed, otherwise the template with the lowest value will be removed.
config.charImage.deblend.maxTempDotProd=0.5

# Apply a smoothing filter to all of the template images
config.charImage.deblend.medianSmoothTemplate=True

# the name of the centroiding algorithm used to set source x,y
config.charImage.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set source moments parameters
config.charImage.measurement.slots.shape='base_SdssShape'

# the name of the algorithm used to set the source aperture flux slot
config.charImage.measurement.slots.apFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source model flux slot
config.charImage.measurement.slots.modelFlux='base_GaussianFlux'

# the name of the algorithm used to set the source psf flux slot
config.charImage.measurement.slots.psfFlux='base_PsfFlux'

# the name of the algorithm used to set the source inst flux slot
config.charImage.measurement.slots.instFlux='base_GaussianFlux'

# the name of the flux measurement algorithm used for calibration
config.charImage.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# When measuring, replace other detected footprints with noise?
config.charImage.measurement.doReplaceWithNoise=True

# How to choose mean and variance of the Gaussian noise we generate?
# Allowed values:
# 	measure	Measure clipped mean and variance from the whole image
# 	meta	Mean = 0, variance = the "BGMEAN" metadata entry
# 	variance	Mean = 0, variance = the image's variance
# 
config.charImage.measurement.noiseReplacer.noiseSource='measure'

# Add ann offset to the generated noise.
config.charImage.measurement.noiseReplacer.noiseOffset=0.0

# The seed multiplier value to use for random number generation
#    >= 1: set the seed deterministically based on exposureId
#       0: fall back to the afw.math.Random default constructor (which uses a seed value of 1)
config.charImage.measurement.noiseReplacer.noiseSeedMultiplier=1

# Prefix to give undeblended plugins
config.charImage.measurement.undeblendedPrefix='undeblended_'

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.charImage.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.charImage.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.charImage.measurement.plugins['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_GaussianCentroid'].doMeasure=True

# Do check that the centroid is contained in footprint.
config.charImage.measurement.plugins['base_GaussianCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.charImage.measurement.plugins['base_GaussianCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.charImage.measurement.plugins['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.charImage.measurement.plugins['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.charImage.measurement.plugins['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.charImage.measurement.plugins['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.charImage.measurement.plugins['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.charImage.measurement.plugins['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.charImage.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.charImage.measurement.plugins['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.charImage.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.charImage.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.charImage.measurement.plugins['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.charImage.measurement.plugins['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.charImage.measurement.plugins['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.charImage.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.charImage.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.charImage.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.charImage.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.charImage.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.charImage.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.charImage.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.charImage.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.charImage.measurement.plugins['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.charImage.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.charImage.measurement.plugins['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.charImage.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.charImage.measurement.plugins['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.charImage.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.charImage.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_SkyCoord'].doMeasure=True

config.charImage.measurement.plugins.names=['base_GaussianFlux', 'base_PixelFlags', 'base_PsfFlux', 'base_SdssCentroid', 'base_SdssShape', 'base_CircularApertureFlux']
# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.charImage.measurement.undeblended['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.charImage.measurement.undeblended['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.charImage.measurement.undeblended['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_GaussianCentroid'].doMeasure=True

# Do check that the centroid is contained in footprint.
config.charImage.measurement.undeblended['base_GaussianCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.charImage.measurement.undeblended['base_GaussianCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.charImage.measurement.undeblended['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.charImage.measurement.undeblended['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.charImage.measurement.undeblended['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.charImage.measurement.undeblended['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.charImage.measurement.undeblended['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.charImage.measurement.undeblended['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.charImage.measurement.undeblended['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.charImage.measurement.undeblended['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.charImage.measurement.undeblended['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.charImage.measurement.undeblended['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.charImage.measurement.undeblended['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.charImage.measurement.undeblended['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.charImage.measurement.undeblended['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.charImage.measurement.undeblended['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.charImage.measurement.undeblended['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.charImage.measurement.undeblended['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.charImage.measurement.undeblended['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.charImage.measurement.undeblended['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.charImage.measurement.undeblended['base_CircularApertureFlux'].maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.charImage.measurement.undeblended['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.charImage.measurement.undeblended['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.charImage.measurement.undeblended['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.charImage.measurement.undeblended['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.charImage.measurement.undeblended['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.charImage.measurement.undeblended['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.charImage.measurement.undeblended['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.charImage.measurement.undeblended['base_Variance'].scale=5.0

# Mask planes to ignore
config.charImage.measurement.undeblended['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_SkyCoord'].doMeasure=True

config.charImage.measurement.undeblended.names=[]
# Run subtasks to measure and apply aperture corrections
config.charImage.doApCorr=True

# Field name prefix for the flux other measurements should be aperture corrected to match
config.charImage.measureApCorr.refFluxName='slot_CalibFlux'

# size of the kernel to create
config.charImage.measureApCorr.starSelector['secondMoment'].kernelSize=21

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measureApCorr.starSelector['secondMoment'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measureApCorr.starSelector['secondMoment'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter']

# specify the minimum psfFlux for good Psf Candidates
config.charImage.measureApCorr.starSelector['secondMoment'].fluxLim=12500.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.charImage.measureApCorr.starSelector['secondMoment'].fluxMax=0.0

# candidate PSF's shapes must lie within this many sigma of the average shape
config.charImage.measureApCorr.starSelector['secondMoment'].clumpNSigma=2.0

# Number of bins in moment histogram
config.charImage.measureApCorr.starSelector['secondMoment'].histSize=64

# Maximum moment to consider
config.charImage.measureApCorr.starSelector['secondMoment'].histMomentMax=100.0

# Multiplier of mean for maximum moments histogram range
config.charImage.measureApCorr.starSelector['secondMoment'].histMomentMaxMultiplier=5.0

# Clipping threshold for moments histogram range
config.charImage.measureApCorr.starSelector['secondMoment'].histMomentClip=5.0

# Multiplier of mean for minimum moments histogram range
config.charImage.measureApCorr.starSelector['secondMoment'].histMomentMinMultiplier=2.0

# size of the kernel to create
config.charImage.measureApCorr.starSelector['objectSize'].kernelSize=21

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measureApCorr.starSelector['objectSize'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measureApCorr.starSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# specify the minimum psfFlux for good Psf Candidates
config.charImage.measureApCorr.starSelector['objectSize'].fluxMin=12500.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.charImage.measureApCorr.starSelector['objectSize'].fluxMax=0.0

# minimum width to include in histogram
config.charImage.measureApCorr.starSelector['objectSize'].widthMin=0.0

# maximum width to include in histogram
config.charImage.measureApCorr.starSelector['objectSize'].widthMax=10.0

# Name of field in Source to use for flux measurement
config.charImage.measureApCorr.starSelector['objectSize'].sourceFluxField='base_GaussianFlux_flux'

# Standard deviation of width allowed to be interpreted as good stars
config.charImage.measureApCorr.starSelector['objectSize'].widthStdAllowed=0.15

# Keep objects within this many sigma of cluster 0's median
config.charImage.measureApCorr.starSelector['objectSize'].nSigmaClip=2.0

# size of the kernel to create
config.charImage.measureApCorr.starSelector['flagged'].kernelSize=21

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measureApCorr.starSelector['flagged'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measureApCorr.starSelector['flagged'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Name of a flag field that is True for stars that should be used.
config.charImage.measureApCorr.starSelector['flagged'].field='calib_psfUsed'

# size of the kernel to create
config.charImage.measureApCorr.starSelector['catalog'].kernelSize=21

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measureApCorr.starSelector['catalog'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measureApCorr.starSelector['catalog'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter']

# specify the minimum psfFlux for good Psf Candidates
# 	Valid Range = [0.0,inf)
config.charImage.measureApCorr.starSelector['catalog'].fluxLim=0.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
# 	Valid Range = [0.0,inf)
config.charImage.measureApCorr.starSelector['catalog'].fluxMax=0.0

config.charImage.measureApCorr.starSelector.name='flagged'
# Minimum number of degrees of freedom (# of valid data points - # of parameters); if this is exceeded, the order of the fit is decreased (in both dimensions), and if we can't decrease it enough, we'll raise ValueError.
# 	Valid Range = [1,inf)
config.charImage.measureApCorr.minDegreesOfFreedom=1

# maximum Chebyshev function order in x
config.charImage.measureApCorr.fitConfig.orderX=2

# maximum Chebyshev function order in y
config.charImage.measureApCorr.fitConfig.orderY=2

# if true, only include terms where the sum of the x and y order is less than or equal to max(orderX, orderY)
config.charImage.measureApCorr.fitConfig.triangular=True

# Number of iterations for sigma clipping
config.charImage.measureApCorr.numIter=4

# Number of standard devisations to clip at
config.charImage.measureApCorr.numSigmaClip=3.0

# Allow these measurement algorithms to fail without an exception
config.charImage.measureApCorr.allowFailure=[]

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.charImage.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.charImage.applyApCorr.doFlagApCorrFailures=True

# flux measurement algorithms to be aperture-corrected by reference to another algorithm; this is a mapping alg1:alg2, where 'alg1' is the algorithm being corrected, and 'alg2' is the algorithm supplying the corrections
config.charImage.applyApCorr.proxies={}

# critical ratio of model to psf flux
config.charImage.catalogCalculation.plugins['base_ClassificationExtendedness'].fluxRatio=0.925

# correction factor for modelFlux error
config.charImage.catalogCalculation.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# correction factor for psfFlux error
config.charImage.catalogCalculation.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

config.charImage.catalogCalculation.plugins.names=['base_ClassificationExtendedness', 'base_FootprintArea']
# Replace the existing PSF model with a simplified version that has the same sigma at the start of each PSF determination iteration? Doing so makes PSF determination converge more robustly and quickly.
config.charImage.useSimplePsf=True

# Estimated FWHM of simple Gaussian PSF model, in pixels. Ignored if input exposure has a PSF model.
config.charImage.installSimplePsf.fwhm=3.5322300675464238

# Width and height of PSF model, in pixels. Must be odd.
# 	Valid Range = [1,inf)
config.charImage.installSimplePsf.width=11

# Padding to add to 4 all edges of the bounding box (pixels)
# 	Valid Range = [0,inf)
config.charImage.refObjLoader.pixelMargin=50

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.charImage.refObjLoader.defaultFilter=''

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.charImage.refObjLoader.filterMap={}

# Maximum separation between reference objects and sources beyond which they will not be considered a match (arcsec)
# 	Valid Range = [0,inf)
config.charImage.ref_match.matcher.maxMatchDistArcSec=3.0

# Number of bright stars to use
# 	Valid Range = [2,inf)
config.charImage.ref_match.matcher.numBrightStars=50

# Minimum number of matched pairs; see also minFracMatchedPairs
# 	Valid Range = [2,inf)
config.charImage.ref_match.matcher.minMatchedPairs=30

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs
# 	Valid Range = [0,1)
config.charImage.ref_match.matcher.minFracMatchedPairs=0.3

# Maximum allowed shift of WCS, due to matching (pixel)
# 	Valid Range = [-inf,4000)
config.charImage.ref_match.matcher.maxOffsetPix=300

# Rotation angle allowed between sources and position reference objects (degrees)
# 	Valid Range = [-inf,6.0)
config.charImage.ref_match.matcher.maxRotationDeg=1.0

# Allowed non-perpendicularity of x and y (degree)
# 	Valid Range = [-inf,45.0)
config.charImage.ref_match.matcher.allowedNonperpDeg=3.0

# number of points to define a shape for matching
config.charImage.ref_match.matcher.numPointsForShape=6

# maximum determinant of linear transformation matrix for a usable solution
config.charImage.ref_match.matcher.maxDeterminant=0.02

# List of flags which cause a source to be rejected as bad
config.charImage.ref_match.matcher.sourceSelector['astrometry'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Type of source flux; typically one of Ap or Psf
config.charImage.ref_match.matcher.sourceSelector['astrometry'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.charImage.ref_match.matcher.sourceSelector['astrometry'].minSnr=10.0

# List of flags which cause a source to be rejected as bad
config.charImage.ref_match.matcher.sourceSelector['matcher'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Type of source flux; typically one of Ap or Psf
config.charImage.ref_match.matcher.sourceSelector['matcher'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.charImage.ref_match.matcher.sourceSelector['matcher'].minSnr=40.0

# List of flags which cause a source to be rejected as bad
config.charImage.ref_match.matcher.sourceSelector['matcherPessimistic'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Type of source flux; typically one of Ap or Psf
config.charImage.ref_match.matcher.sourceSelector['matcherPessimistic'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.charImage.ref_match.matcher.sourceSelector['matcherPessimistic'].minSnr=40.0

config.charImage.ref_match.matcher.sourceSelector.name='matcher'
# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.charImage.ref_match.matchDistanceSigma=2.0

# size of the kernel to create
config.charImage.measurePsf.starSelector['secondMoment'].kernelSize=21

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measurePsf.starSelector['secondMoment'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measurePsf.starSelector['secondMoment'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter']

# specify the minimum psfFlux for good Psf Candidates
config.charImage.measurePsf.starSelector['secondMoment'].fluxLim=12500.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.charImage.measurePsf.starSelector['secondMoment'].fluxMax=0.0

# candidate PSF's shapes must lie within this many sigma of the average shape
config.charImage.measurePsf.starSelector['secondMoment'].clumpNSigma=2.0

# Number of bins in moment histogram
config.charImage.measurePsf.starSelector['secondMoment'].histSize=64

# Maximum moment to consider
config.charImage.measurePsf.starSelector['secondMoment'].histMomentMax=100.0

# Multiplier of mean for maximum moments histogram range
config.charImage.measurePsf.starSelector['secondMoment'].histMomentMaxMultiplier=5.0

# Clipping threshold for moments histogram range
config.charImage.measurePsf.starSelector['secondMoment'].histMomentClip=5.0

# Multiplier of mean for minimum moments histogram range
config.charImage.measurePsf.starSelector['secondMoment'].histMomentMinMultiplier=2.0

# size of the kernel to create
config.charImage.measurePsf.starSelector['objectSize'].kernelSize=21

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measurePsf.starSelector['objectSize'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measurePsf.starSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# specify the minimum psfFlux for good Psf Candidates
config.charImage.measurePsf.starSelector['objectSize'].fluxMin=12500.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.charImage.measurePsf.starSelector['objectSize'].fluxMax=0.0

# minimum width to include in histogram
config.charImage.measurePsf.starSelector['objectSize'].widthMin=0.0

# maximum width to include in histogram
config.charImage.measurePsf.starSelector['objectSize'].widthMax=10.0

# Name of field in Source to use for flux measurement
config.charImage.measurePsf.starSelector['objectSize'].sourceFluxField='base_GaussianFlux_flux'

# Standard deviation of width allowed to be interpreted as good stars
config.charImage.measurePsf.starSelector['objectSize'].widthStdAllowed=0.15

# Keep objects within this many sigma of cluster 0's median
config.charImage.measurePsf.starSelector['objectSize'].nSigmaClip=2.0

# size of the kernel to create
config.charImage.measurePsf.starSelector['flagged'].kernelSize=21

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measurePsf.starSelector['flagged'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measurePsf.starSelector['flagged'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Name of a flag field that is True for stars that should be used.
config.charImage.measurePsf.starSelector['flagged'].field='calib_psfUsed'

# size of the kernel to create
config.charImage.measurePsf.starSelector['catalog'].kernelSize=21

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measurePsf.starSelector['catalog'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measurePsf.starSelector['catalog'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter']

# specify the minimum psfFlux for good Psf Candidates
# 	Valid Range = [0.0,inf)
config.charImage.measurePsf.starSelector['catalog'].fluxLim=0.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
# 	Valid Range = [0.0,inf)
config.charImage.measurePsf.starSelector['catalog'].fluxMax=0.0

config.charImage.measurePsf.starSelector.name='objectSize'
# radius of the kernel to create, relative to the square root of the stellar quadrupole moments
config.charImage.measurePsf.psfDeterminer['pca'].kernelSize=10.0

# Minimum radius of the kernel
config.charImage.measurePsf.psfDeterminer['pca'].kernelSizeMin=25

# Maximum radius of the kernel
config.charImage.measurePsf.psfDeterminer['pca'].kernelSizeMax=45

# Use non-linear fitter for spatial variation of Kernel
config.charImage.measurePsf.psfDeterminer['pca'].nonLinearSpatialFit=False

# number of eigen components for PSF kernel creation
config.charImage.measurePsf.psfDeterminer['pca'].nEigenComponents=4

# specify spatial order for PSF kernel creation
config.charImage.measurePsf.psfDeterminer['pca'].spatialOrder=2

# size of cell used to determine PSF (pixels, column direction)
config.charImage.measurePsf.psfDeterminer['pca'].sizeCellX=256

# size of cell used to determine PSF (pixels, row direction)
config.charImage.measurePsf.psfDeterminer['pca'].sizeCellY=256

# number of stars per psf cell for PSF kernel creation
config.charImage.measurePsf.psfDeterminer['pca'].nStarPerCell=3

# Number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measurePsf.psfDeterminer['pca'].borderWidth=0

# number of stars per psf Cell for spatial fitting
config.charImage.measurePsf.psfDeterminer['pca'].nStarPerCellSpatialFit=5

# Should each PSF candidate be given the same weight, independent of magnitude?
config.charImage.measurePsf.psfDeterminer['pca'].constantWeight=True

# number of iterations of PSF candidate star list
config.charImage.measurePsf.psfDeterminer['pca'].nIterForPsf=3

# tolerance of spatial fitting
config.charImage.measurePsf.psfDeterminer['pca'].tolerance=0.01

# floor for variance is lam*data
config.charImage.measurePsf.psfDeterminer['pca'].lam=0.05

# for psf candidate evaluation
config.charImage.measurePsf.psfDeterminer['pca'].reducedChi2ForPsfCandidates=2.0

# Rejection threshold (stdev) for candidates based on spatial fit
config.charImage.measurePsf.psfDeterminer['pca'].spatialReject=3.0

# Threshold (stdev) for rejecting extraneous pixels around candidate; applied if positive
config.charImage.measurePsf.psfDeterminer['pca'].pixelThreshold=0.0

# Reject candidates that are blended?
config.charImage.measurePsf.psfDeterminer['pca'].doRejectBlends=False

# Mask blends in image?
config.charImage.measurePsf.psfDeterminer['pca'].doMaskBlends=True

config.charImage.measurePsf.psfDeterminer.name='pca'
# Fraction of PSF candidates to reserve from fitting; none if <= 0
config.charImage.measurePsf.reserveFraction=-1.0

# This number will be multiplied by the exposure ID to set the random seed for reserving candidates
config.charImage.measurePsf.reserveSeed=1

# Interpolate over defects? (ignored unless you provide a list of defects)
config.charImage.repair.doInterpolate=True

# Find and mask out cosmic rays?
config.charImage.repair.doCosmicRay=True

# maximum number of contaminated pixels
config.charImage.repair.cosmicray.nCrPixelMax=100000

# CRs must be > this many sky-sig above sky
config.charImage.repair.cosmicray.minSigma=6.0

# CRs must have > this many DN (== electrons/gain) in initial detection
config.charImage.repair.cosmicray.min_DN=150.0

# used in condition 3 for CR; see CR.cc code
config.charImage.repair.cosmicray.cond3_fac=2.5

# used in condition 3 for CR; see CR.cc code
config.charImage.repair.cosmicray.cond3_fac2=0.6

# number of times to look for contaminated pixels near known CR pixels
config.charImage.repair.cosmicray.niteration=3

# Don't interpolate over CR pixels
config.charImage.repair.cosmicray.keepCRs=False

# type of statistic to use for grid points
# Allowed values:
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 	MEDIAN	median
# 	None	Field is optional
# 
config.charImage.repair.cosmicray.background.statisticsProperty='MEDIAN'

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 
config.charImage.repair.cosmicray.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.charImage.repair.cosmicray.background.binSize=100000

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
# 	Valid Range = [0,inf)
config.charImage.repair.cosmicray.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
# 	Valid Range = [0,inf)
config.charImage.repair.cosmicray.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NONE	No background estimation is to be attempted
# 	None	Field is optional
# 
config.charImage.repair.cosmicray.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.charImage.repair.cosmicray.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.charImage.repair.cosmicray.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.charImage.repair.cosmicray.background.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.charImage.repair.cosmicray.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.charImage.repair.cosmicray.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.charImage.repair.cosmicray.background.weighting=True

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.charImage.repair.interp.modelPsf.size=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.charImage.repair.interp.modelPsf.sizeFactor=3.0

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.charImage.repair.interp.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.charImage.repair.interp.modelPsf.maxSize=None

# Default FWHM of Gaussian model of core of star (pixels)
config.charImage.repair.interp.modelPsf.defaultFwhm=3.0

# Add a Gaussian to represent wings?
config.charImage.repair.interp.modelPsf.addWing=True

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.charImage.repair.interp.modelPsf.wingFwhmFactor=2.5

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.charImage.repair.interp.modelPsf.wingAmplitude=0.1

# Smoothly taper to the fallback value at the edge of the image?
config.charImage.repair.interp.useFallbackValueAtEdge=True

# Type of statistic to calculate edge fallbackValue for interpolation
# Allowed values:
# 	MEAN	mean
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	USER	user value set in fallbackUserValue config
# 	None	Field is optional
# 
config.charImage.repair.interp.fallbackValueType='MEANCLIP'

# If fallbackValueType is 'USER' then use this as the fallbackValue; ignored otherwise
config.charImage.repair.interp.fallbackUserValue=0.0

# Allow negative values for egde interpolation fallbackValue?  If False, set fallbackValue to max(fallbackValue, 0.0)
config.charImage.repair.interp.negativeFallbackAllowed=True

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.charImage.checkUnitsParseStrict='raise'

# Perform calibration?
config.doCalibrate=True

# Save calibration results?
config.calibrate.doWrite=True

# Include HeavyFootprint data in source table? If false then heavy footprints are saved as normal footprints, which saves some space
config.calibrate.doWriteHeavyFootprintsInSources=True

# Write reference matches (ignored if doWrite false)?
config.calibrate.doWriteMatches=True

# Write reference matches in denormalized format? This format uses more disk space, but is more convenient to read. Ignored if doWriteMatches=False or doWrite=False.
config.calibrate.doWriteMatchesDenormalized=False

# Perform astrometric calibration?
config.calibrate.doAstrometry=True

# Padding to add to 4 all edges of the bounding box (pixels)
# 	Valid Range = [0,inf)
config.calibrate.astromRefObjLoader.pixelMargin=50

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.calibrate.astromRefObjLoader.defaultFilter=''

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.calibrate.astromRefObjLoader.filterMap={}

# Padding to add to 4 all edges of the bounding box (pixels)
# 	Valid Range = [0,inf)
config.calibrate.photoRefObjLoader.pixelMargin=50

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.calibrate.photoRefObjLoader.defaultFilter=''

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.calibrate.photoRefObjLoader.filterMap={}

# Maximum separation between reference objects and sources beyond which they will not be considered a match (arcsec)
# 	Valid Range = [0,inf)
config.calibrate.astrometry.matcher.maxMatchDistArcSec=3.0

# Number of bright stars to use
# 	Valid Range = [2,inf)
config.calibrate.astrometry.matcher.numBrightStars=50

# Minimum number of matched pairs; see also minFracMatchedPairs
# 	Valid Range = [2,inf)
config.calibrate.astrometry.matcher.minMatchedPairs=30

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs
# 	Valid Range = [0,1)
config.calibrate.astrometry.matcher.minFracMatchedPairs=0.3

# Maximum allowed shift of WCS, due to matching (pixel)
# 	Valid Range = [-inf,4000)
config.calibrate.astrometry.matcher.maxOffsetPix=300

# Rotation angle allowed between sources and position reference objects (degrees)
# 	Valid Range = [-inf,6.0)
config.calibrate.astrometry.matcher.maxRotationDeg=1.0

# Allowed non-perpendicularity of x and y (degree)
# 	Valid Range = [-inf,45.0)
config.calibrate.astrometry.matcher.allowedNonperpDeg=3.0

# number of points to define a shape for matching
config.calibrate.astrometry.matcher.numPointsForShape=6

# maximum determinant of linear transformation matrix for a usable solution
config.calibrate.astrometry.matcher.maxDeterminant=0.02

# List of flags which cause a source to be rejected as bad
config.calibrate.astrometry.matcher.sourceSelector['astrometry'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Type of source flux; typically one of Ap or Psf
config.calibrate.astrometry.matcher.sourceSelector['astrometry'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.calibrate.astrometry.matcher.sourceSelector['astrometry'].minSnr=10.0

# List of flags which cause a source to be rejected as bad
config.calibrate.astrometry.matcher.sourceSelector['matcher'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Type of source flux; typically one of Ap or Psf
config.calibrate.astrometry.matcher.sourceSelector['matcher'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.calibrate.astrometry.matcher.sourceSelector['matcher'].minSnr=40.0

# List of flags which cause a source to be rejected as bad
config.calibrate.astrometry.matcher.sourceSelector['matcherPessimistic'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Type of source flux; typically one of Ap or Psf
config.calibrate.astrometry.matcher.sourceSelector['matcherPessimistic'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.calibrate.astrometry.matcher.sourceSelector['matcherPessimistic'].minSnr=40.0

config.calibrate.astrometry.matcher.sourceSelector.name='matcher'
# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.calibrate.astrometry.matchDistanceSigma=2.0

# order of SIP polynomial
# 	Valid Range = [0,inf)
config.calibrate.astrometry.wcsFitter.order=4

# number of iterations of fitter (which fits X and Y separately, and so benefits from a few iterations
# 	Valid Range = [1,inf)
config.calibrate.astrometry.wcsFitter.numIter=3

# number of rejection iterations
# 	Valid Range = [0,inf)
config.calibrate.astrometry.wcsFitter.numRejIter=1

# Number of standard deviations for clipping level
# 	Valid Range = [0.0,inf)
config.calibrate.astrometry.wcsFitter.rejSigma=3.0

# maximum median scatter of a WCS fit beyond which the fit fails (arcsec); be generous, as this is only intended to catch catastrophic failures
# 	Valid Range = [0,inf)
config.calibrate.astrometry.wcsFitter.maxScatterArcsec=10.0

# If True then load reference objects and match sources but do not fit a WCS;  this simply controls whether 'run' calls 'solve' or 'loadAndMatch'
config.calibrate.astrometry.forceKnownWcs=False

# maximum number of iterations of match sources and fit WCSignored if not fitting a WCS
# 	Valid Range = [1,inf)
config.calibrate.astrometry.maxIter=3

# the match distance below which further iteration is pointless (arcsec); ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.calibrate.astrometry.minMatchDistanceArcSec=0.001

# Raise an exception if astrometry fails? Ignored if doAstrometry false.
config.calibrate.requireAstrometry=True

# Perform phometric calibration?
config.calibrate.doPhotoCal=True

# Raise an exception if photoCal fails? Ignored if doPhotoCal false.
config.calibrate.requirePhotoCal=True

# Maximum separation between reference objects and sources beyond which they will not be considered a match (arcsec)
# 	Valid Range = [0,inf)
config.calibrate.photoCal.matcher.maxMatchDistArcSec=3.0

# Number of bright stars to use
# 	Valid Range = [2,inf)
config.calibrate.photoCal.matcher.numBrightStars=50

# Minimum number of matched pairs; see also minFracMatchedPairs
# 	Valid Range = [2,inf)
config.calibrate.photoCal.matcher.minMatchedPairs=30

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs
# 	Valid Range = [0,1)
config.calibrate.photoCal.matcher.minFracMatchedPairs=0.3

# Maximum allowed shift of WCS, due to matching (pixel)
# 	Valid Range = [-inf,4000)
config.calibrate.photoCal.matcher.maxOffsetPix=300

# Rotation angle allowed between sources and position reference objects (degrees)
# 	Valid Range = [-inf,6.0)
config.calibrate.photoCal.matcher.maxRotationDeg=1.0

# Allowed non-perpendicularity of x and y (degree)
# 	Valid Range = [-inf,45.0)
config.calibrate.photoCal.matcher.allowedNonperpDeg=3.0

# number of points to define a shape for matching
config.calibrate.photoCal.matcher.numPointsForShape=6

# maximum determinant of linear transformation matrix for a usable solution
config.calibrate.photoCal.matcher.maxDeterminant=0.02

# List of flags which cause a source to be rejected as bad
config.calibrate.photoCal.matcher.sourceSelector['astrometry'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Type of source flux; typically one of Ap or Psf
config.calibrate.photoCal.matcher.sourceSelector['astrometry'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.calibrate.photoCal.matcher.sourceSelector['astrometry'].minSnr=10.0

# List of flags which cause a source to be rejected as bad
config.calibrate.photoCal.matcher.sourceSelector['matcher'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Type of source flux; typically one of Ap or Psf
config.calibrate.photoCal.matcher.sourceSelector['matcher'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.calibrate.photoCal.matcher.sourceSelector['matcher'].minSnr=40.0

# List of flags which cause a source to be rejected as bad
config.calibrate.photoCal.matcher.sourceSelector['matcherPessimistic'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Type of source flux; typically one of Ap or Psf
config.calibrate.photoCal.matcher.sourceSelector['matcherPessimistic'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.calibrate.photoCal.matcher.sourceSelector['matcherPessimistic'].minSnr=40.0

config.calibrate.photoCal.matcher.sourceSelector.name='matcher'
# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.calibrate.photoCal.matchDistanceSigma=2.0

# Don't use objects fainter than this magnitude
config.calibrate.photoCal.magLimit=22.0

# Fraction of candidates to reserve from fitting; none if <= 0
config.calibrate.photoCal.reserveFraction=-1.0

# This number will be multiplied by the exposure ID to set the random seed for reserving candidates
config.calibrate.photoCal.reserveSeed=1

# Name of the source flux field to use.  The associated flag field
# ('<name>_flags') will be implicitly included in badFlags.
config.calibrate.photoCal.fluxField='slot_CalibFlux_flux'

# Apply photometric color terms to reference stars? One of:
# None: apply if colorterms and photoCatName are not None;
#       fail if color term data is not available for the specified ref catalog and filter.
# True: always apply colorterms; fail if color term data is not available for the
#       specified reference catalog and filter.
# False: do not apply.
config.calibrate.photoCal.applyColorTerms=None

# List of source flag fields that must be set for a source to be used.
config.calibrate.photoCal.goodFlags=[]

# List of source flag fields that will cause a source to be rejected when they are set.
config.calibrate.photoCal.badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolated', 'base_PixelFlags_flag_saturated']

# maximum sigma to use when clipping
config.calibrate.photoCal.sigmaMax=0.25

# clip at nSigma
config.calibrate.photoCal.nSigma=3.0

# use median instead of mean to compute zeropoint
config.calibrate.photoCal.useMedian=True

# number of iterations
config.calibrate.photoCal.nIter=20

config.calibrate.photoCal.colorterms.data={}
# Name of photometric reference catalog; used to select a color term dict in colorterms. see also applyColorTerms
config.calibrate.photoCal.photoCatName='10.0+111'

# Additional magnitude uncertainty to be added in quadrature with measurement errors.
# 	Valid Range = [0.0,inf)
config.calibrate.photoCal.magErrFloor=0.0

# Use the extendedness parameter to select objects to use in photometric calibration?
# This applies only to the sources detected on the exposure, not the reference catalog
config.calibrate.photoCal.doSelectUnresolved=True

# Fields to copy from the icSource catalog to the output catalog for matching sources Any missing fields will trigger a RuntimeError exception. Ignored if icSourceCat is not provided.
config.calibrate.icSourceFieldsToCopy=['calib_psfCandidate', 'calib_psfUsed', 'calib_psfReserved']

# Match radius for matching icSourceCat objects to sourceCat objects (pixels)
config.calibrate.matchRadiusPix=3.0

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.calibrate.checkUnitsParseStrict='raise'

# detected sources with fewer than the specified number of pixels will be ignored
# 	Valid Range = [0,inf)
config.calibrate.detection.minPixels=1

# Pixels should be grown as isotropically as possible (slower)
config.calibrate.detection.isotropicGrow=False

# Grow detections by nSigmaToGrow * sigma; if 0 then do not grow
config.calibrate.detection.nSigmaToGrow=2.4

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.calibrate.detection.returnOriginalFootprints=False

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.calibrate.detection.thresholdValue=5.0

# Include threshold relative to thresholdValue
# 	Valid Range = [0.0,inf)
config.calibrate.detection.includeThresholdMultiplier=1.0

# specifies the desired flavor of Threshold
# Allowed values:
# 	variance	threshold applied to image variance
# 	stdev	threshold applied to image std deviation
# 	value	threshold applied to image value
# 	pixel_stdev	threshold applied to per-pixel std deviation
# 
config.calibrate.detection.thresholdType='stdev'

# specifies whether to detect positive, or negative sources, or both
# Allowed values:
# 	positive	detect only positive sources
# 	negative	detect only negative sources
# 	both	detect both positive and negative sources
# 
config.calibrate.detection.thresholdPolarity='positive'

# Fiddle factor to add to the background; debugging only
config.calibrate.detection.adjustBackground=0.0

# Estimate the background again after final source detection?
config.calibrate.detection.reEstimateBackground=True

# type of statistic to use for grid points
# Allowed values:
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 	MEDIAN	median
# 	None	Field is optional
# 
config.calibrate.detection.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 
config.calibrate.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.calibrate.detection.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
# 	Valid Range = [0,inf)
config.calibrate.detection.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
# 	Valid Range = [0,inf)
config.calibrate.detection.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NONE	No background estimation is to be attempted
# 	None	Field is optional
# 
config.calibrate.detection.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.calibrate.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.calibrate.detection.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.calibrate.detection.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.calibrate.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.calibrate.detection.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.calibrate.detection.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 	MEDIAN	median
# 	None	Field is optional
# 
config.calibrate.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 
config.calibrate.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.calibrate.detection.tempLocalBackground.binSize=64

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
# 	Valid Range = [0,inf)
config.calibrate.detection.tempLocalBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
# 	Valid Range = [0,inf)
config.calibrate.detection.tempLocalBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NONE	No background estimation is to be attempted
# 	None	Field is optional
# 
config.calibrate.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.calibrate.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.calibrate.detection.tempLocalBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.calibrate.detection.tempLocalBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.calibrate.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.calibrate.detection.tempLocalBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.calibrate.detection.tempLocalBackground.weighting=True

# Do temporary interpolated background subtraction before footprint detection?
config.calibrate.detection.doTempLocalBackground=False

# The maximum number of peaks in a Footprint before trying to replace its peaks using the temporary local background
config.calibrate.detection.nPeaksMaxSimple=1

# Run deblender input exposure
config.calibrate.doDeblend=True

# What to do when a peak to be deblended is close to the edge of the image
# Allowed values:
# 	clip	Clip the template at the edge AND the mirror of the edge.
# 	ramp	Ramp down flux at the image edge by the PSF
# 	noclip	Ignore the edge when building the symmetric template.
# 	None	Field is optional
# 
config.calibrate.deblend.edgeHandling='ramp'

# When the deblender should attribute stray flux to point sources
# Allowed values:
# 	necessary	When there is not an extended object in the footprint
# 	always	Always
# 	never	Never; stray flux will not be attributed to any deblended child if the deblender thinks all peaks look like point sources
# 	None	Field is optional
# 
config.calibrate.deblend.strayFluxToPointSources='necessary'

# Assign stray flux (not claimed by any child in the deblender) to deblend children.
config.calibrate.deblend.assignStrayFlux=True

# How to split flux among peaks
# Allowed values:
# 	r-to-peak	~ 1/(1+R^2) to the peak
# 	r-to-footprint	~ 1/(1+R^2) to the closest pixel in the footprint.  CAUTION: this can be computationally expensive on large footprints!
# 	nearest-footprint	Assign 100% to the nearest footprint (using L-1 norm aka Manhattan distance)
# 	trim	Shrink the parent footprint to pixels that are not assigned to children
# 	None	Field is optional
# 
config.calibrate.deblend.strayFluxRule='trim'

# When splitting stray flux, clip fractions below this value to zero.
config.calibrate.deblend.clipStrayFluxFraction=0.001

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.calibrate.deblend.psfChisq1=1.5

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.calibrate.deblend.psfChisq2=1.5

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.calibrate.deblend.psfChisq2b=1.5

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.calibrate.deblend.maxNumberOfPeaks=0

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.calibrate.deblend.maxFootprintArea=1000000

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.calibrate.deblend.maxFootprintSize=0

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.calibrate.deblend.minFootprintAxisRatio=0.0

# Mask name for footprints not deblended, or None
config.calibrate.deblend.notDeblendedMask='NOT_DEBLENDED'

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
# 	Valid Range = [2,inf)
config.calibrate.deblend.tinyFootprintSize=2

# Guarantee that all peaks produce a child source.
config.calibrate.deblend.propagateAllPeaks=False

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.calibrate.deblend.catchFailures=False

# Mask planes to ignore when performing statistics
config.calibrate.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.calibrate.deblend.maskLimits={}

# If true, a least-squares fit of the templates will be done to the full image. The templates will be re-weighted based on this fit.
config.calibrate.deblend.weightTemplates=False

# Try to remove similar templates?
config.calibrate.deblend.removeDegenerateTemplates=False

# If the dot product between two templates is larger than this value, we consider them to be describing the same object (i.e. they are degenerate).  If one of the objects has been labeled as a PSF it will be removed, otherwise the template with the lowest value will be removed.
config.calibrate.deblend.maxTempDotProd=0.5

# Apply a smoothing filter to all of the template images
config.calibrate.deblend.medianSmoothTemplate=True

# the name of the centroiding algorithm used to set source x,y
config.calibrate.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set source moments parameters
config.calibrate.measurement.slots.shape='base_SdssShape'

# the name of the algorithm used to set the source aperture flux slot
config.calibrate.measurement.slots.apFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source model flux slot
config.calibrate.measurement.slots.modelFlux='base_GaussianFlux'

# the name of the algorithm used to set the source psf flux slot
config.calibrate.measurement.slots.psfFlux='base_PsfFlux'

# the name of the algorithm used to set the source inst flux slot
config.calibrate.measurement.slots.instFlux='base_GaussianFlux'

# the name of the flux measurement algorithm used for calibration
config.calibrate.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# When measuring, replace other detected footprints with noise?
config.calibrate.measurement.doReplaceWithNoise=True

# How to choose mean and variance of the Gaussian noise we generate?
# Allowed values:
# 	measure	Measure clipped mean and variance from the whole image
# 	meta	Mean = 0, variance = the "BGMEAN" metadata entry
# 	variance	Mean = 0, variance = the image's variance
# 
config.calibrate.measurement.noiseReplacer.noiseSource='measure'

# Add ann offset to the generated noise.
config.calibrate.measurement.noiseReplacer.noiseOffset=0.0

# The seed multiplier value to use for random number generation
#    >= 1: set the seed deterministically based on exposureId
#       0: fall back to the afw.math.Random default constructor (which uses a seed value of 1)
config.calibrate.measurement.noiseReplacer.noiseSeedMultiplier=1

# Prefix to give undeblended plugins
config.calibrate.measurement.undeblendedPrefix='undeblended_'

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.calibrate.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.calibrate.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.calibrate.measurement.plugins['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_GaussianCentroid'].doMeasure=True

# Do check that the centroid is contained in footprint.
config.calibrate.measurement.plugins['base_GaussianCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.calibrate.measurement.plugins['base_GaussianCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.calibrate.measurement.plugins['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.calibrate.measurement.plugins['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.calibrate.measurement.plugins['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.calibrate.measurement.plugins['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.calibrate.measurement.plugins['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.calibrate.measurement.plugins['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.calibrate.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.calibrate.measurement.plugins['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.calibrate.measurement.plugins['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.calibrate.measurement.plugins['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.calibrate.measurement.plugins['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.calibrate.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.calibrate.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.calibrate.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.calibrate.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.calibrate.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.calibrate.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.calibrate.measurement.plugins['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.calibrate.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.calibrate.measurement.plugins['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.calibrate.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.calibrate.measurement.plugins['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.calibrate.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.calibrate.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_SkyCoord'].doMeasure=True

config.calibrate.measurement.plugins.names=['base_GaussianFlux', 'base_PixelFlags', 'base_PsfFlux', 'base_GaussianCentroid', 'base_Blendedness', 'base_NaiveCentroid', 'base_Variance', 'base_SkyCoord', 'base_SdssCentroid', 'base_SdssShape', 'base_CircularApertureFlux']
# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.calibrate.measurement.undeblended['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.calibrate.measurement.undeblended['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.calibrate.measurement.undeblended['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_GaussianCentroid'].doMeasure=True

# Do check that the centroid is contained in footprint.
config.calibrate.measurement.undeblended['base_GaussianCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.calibrate.measurement.undeblended['base_GaussianCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.calibrate.measurement.undeblended['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.calibrate.measurement.undeblended['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.calibrate.measurement.undeblended['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.calibrate.measurement.undeblended['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.calibrate.measurement.undeblended['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.calibrate.measurement.undeblended['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.calibrate.measurement.undeblended['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.calibrate.measurement.undeblended['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.measurement.undeblended['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.measurement.undeblended['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.calibrate.measurement.undeblended['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.calibrate.measurement.undeblended['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.calibrate.measurement.undeblended['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.calibrate.measurement.undeblended['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.calibrate.measurement.undeblended['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.calibrate.measurement.undeblended['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.calibrate.measurement.undeblended['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.measurement.undeblended['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.calibrate.measurement.undeblended['base_CircularApertureFlux'].maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.calibrate.measurement.undeblended['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.measurement.undeblended['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.calibrate.measurement.undeblended['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.calibrate.measurement.undeblended['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.calibrate.measurement.undeblended['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.calibrate.measurement.undeblended['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.calibrate.measurement.undeblended['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.calibrate.measurement.undeblended['base_Variance'].scale=5.0

# Mask planes to ignore
config.calibrate.measurement.undeblended['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_SkyCoord'].doMeasure=True

config.calibrate.measurement.undeblended.names=[]
# Run subtask to apply aperture correction
config.calibrate.doApCorr=True

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.calibrate.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.calibrate.applyApCorr.doFlagApCorrFailures=True

# flux measurement algorithms to be aperture-corrected by reference to another algorithm; this is a mapping alg1:alg2, where 'alg1' is the algorithm being corrected, and 'alg2' is the algorithm supplying the corrections
config.calibrate.applyApCorr.proxies={}

# critical ratio of model to psf flux
config.calibrate.catalogCalculation.plugins['base_ClassificationExtendedness'].fluxRatio=0.925

# correction factor for modelFlux error
config.calibrate.catalogCalculation.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# correction factor for psfFlux error
config.calibrate.catalogCalculation.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

config.calibrate.catalogCalculation.plugins.names=['base_ClassificationExtendedness', 'base_FootprintArea']
# Run fake sources injection task
config.calibrate.doInsertFakes=False

# Mask plane to set on pixels affected by fakes.  Will be added if not already present.
config.calibrate.insertFakes.maskPlaneName='FAKE'

