# use "instcal" exposures from the community pipeline until DecamIsrTask is up to snuff
from lsst.obs.decam.decamNullIsr import DecamNullIsrTask
config.isr.retarget(DecamNullIsrTask)

from lsst.meas.astrom.matchPessimisticB import MatchPessimisticBTask
config.calibrate.astrometry.matcher.retarget(MatchPessimisticBTask)
