## @knitr correct-pitch

XformLB <- function (bvector, .roll, .pitch, .heading, .reverse=FALSE) {
  Cradeg <- pi/180
  sr <- sin(.roll*Cradeg); cr <- cos(.roll*Cradeg)
  sp <- sin(.pitch*Cradeg); cp <- cos(.pitch*Cradeg)
  sh <- sin(.heading*Cradeg); ch <- cos(.heading*Cradeg)
  M <- c (ch*cr+sh*sr*sp, -sh*cr+ch*sr*sp, -cp*sr, sh*cp, ch*cp, sp, 
          ch*sr-sh*sp*cr, -sh*sr-ch*sp*sr, cp*cr)
  dim (M) <- c(3,3)
  if (.reverse) {M <- t(M)}
  return (M %*% bvector)
}

## The correction should be subtracted from PITCH to get PITCHC;
## it is the error in pitch, so you get the true value by subtraction
## This calculates corrections for an entire flight in one call.
## D must be a dataframe containing at least VNS, VEW, GGVNS, GGVEW,
## LATC, GGALT, THDG, PITCH, ROLL, 
CorrectPitch <- function (D) {
  Cradeg = pi/180
  .vns <- zoo::na.approx (as.vector(D$VNS), maxgap=1000, na.rm = FALSE)
  .vew <- zoo::na.approx (as.vector(D$VEW), maxgap=1000, na.rm = FALSE)
  .ggvns <- zoo::na.approx (as.vector(D$GGVNS), maxgap=1000, na.rm = FALSE)
  .ggvew <- zoo::na.approx (as.vector(D$GGVEW), maxgap=1000, na.rm = FALSE)
  # 1013 points (must be odd) to span about 1/5 Schuler osc. -- about 16.8 min
  vndot <- signal::sgolayfilt (.vns-.ggvns, 3, 1013, m=1)  # m=1 for first deriv.
  vedot <- signal::sgolayfilt (.vew-.ggvew, 3, 1013, m=1)
  deltaPitchL <- -vndot/Ranadu::Gravity (D$LATC, D$GGALT)
  deltaRollL  <- -vedot/Ranadu::Gravity (D$LATC, D$GGALT)
  LHDG <- length(D$THDG)
  # .hdg <- D$THDG*Cradeg
  # deltaPitch <- (sin(.hdg)*deltaRollL + cos(.hdg)*deltaPitchL)/Cradeg
  # deltaRoll <- (cos(.hdg)*deltaRollL - sin(.hdg)*deltaPitchL)/Cradeg
  ## replace with the full transformation matrix:
  deltaPitch2 <- vector ("numeric", LHDG)
  deltaRoll2  <- vector ("numeric", LHDG)
  for (i in 1:LHDG) {
    bl <- c(sin(deltaRollL[i]), sin(deltaPitchL[i]), 
            sqrt(1.-sin(deltaRollL[i])^2 - sin(deltaPitchL[i])^2))
    bb <- as.vector(XformLB (bl, D$ROLL[i], D$PITCH[i], D$THDG[i], .reverse=TRUE))
    deltaPitch2[i] <- atan(bb[2]/bb[3]) / Cradeg - D$PITCH[i]
    deltaRoll2[i]  <- atan(bb[1]/bb[3]) / Cradeg - D$ROLL[i]
  }
  #return (deltaPitch)
  return (deltaPitch2)
}

