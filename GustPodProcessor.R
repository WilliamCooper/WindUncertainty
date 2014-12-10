
## ----Initialization,echo=F,include=F-------------------------------------
#' @title GustPodProcessor
#' @description Documents the processing for the gust pod and also implements that processing to add gust-pod variables WDG, WSG, WIG, TASG (analogous to WDC, WSC, WIC, TASX) for a specified project and flight.
#' @details Constructs a memo with the math basis for the transformations required to calculate wind from the gust pod, incorporating calibration for that sensor, and also processes a specified netCDF file to add new variables to that file. The processing uses the angle of attack and sideslip angle as determined from the gust pod and also the true airspeed from that sensor, so produces a measurement that is independent of the primary radome-based system. The memo constructed by running this script also provides plots and other characteristics of the new variables. 
#' @author Al Cooper
#' @export 
#' @param 
#' @return
require(Ranadu)
library(knitr)
Flight <- "rf16"
Directory <- "/home/Data/"
Project <- "DEEPWAVE"
fname = sprintf("%s%s/DEEPWAVE%s.nc", Directory, Project, Flight)
# copy to a new file before adding variables to it:
fnew = sprintf("%s%s/DW%s_GP.nc", Directory, Project,Flight)
system(paste ("cp", fname, fnew, sep=' '), wait=TRUE)
VarNames <- c("VYC","GGALT","LATC", "LONC", "PSXC", "QCXC",
              "WDC", "WSC", "GGVEW", "GGVNS", "VEW", "VNS",
              "ADIFR", "AKRD", "SSLIP", "PITCH", "TASX",
              "ROLL", "THDG", "BDIFR", "EWX",
              "ADIFR_GP", "BDIFR_GP", "PS_GP", "QC_GP",
              "CROLL_GP", "CPITCH_GP", "CTHDG_GP", "WIC",
              "CVNS_GP", "CVEW_GP", "VSPD", "CVSPD_GP",
              "ATX")
Data <- getNetCDF (fname, VarNames) # high-rate data OK here
Cradeg <- pi / 180.



## ----Relative-Wind-Code,echo=F,include=F---------------------------------


attach(Data)    # get the unrestricted full data file
AQR <- ADIFR / QCXC
RR <- QCXC/PSXC
Mach <- MachNumber (PSXC, QCXC)
AQR_GP <- ADIFR_GP/QC_GP
RR_GP <- QC_GP/PS_GP
# use coefficients saved from CalibrationPart1.rnw:
cf_GP <- c(-0.9151291, 3.9277735, 3.1905422, 1.3502541)
cfr <- c(4.402532, 21.872829)
AOA <- cfr[1] + AQR * cfr[2] 
AOA_GP <- cf_GP[1] + AQR_GP * (cf_GP[2] + cf_GP[3]*Mach) + cf_GP[4]*RR_GP
cfs <- c(0., 21.335) 
cfs_GP <- c(-3.498, 11.443)
SS <- cfs[1] + cfs[2] * BDIFR / QCXC
SS_GP <- cfs_GP[1] + cfs_GP[2] * BDIFR_GP / QC_GP
bf <- c(-0.07791, -1.1571, 2.6691, -1.3588, -0.0046686, 0.53110)
RRS_GP <- QC_GP / PS_GP
RRS <- bf[1] + RRS_GP * (bf[2] + bf[3] * Mach + bf[4] * RRS_GP) + bf[5] * ADIFR_GP / QC_GP + bf[6] * Mach
# note that PSXC in the next equation doesn't matter, only the ratio RRS
TAS_GP <- TrueAirspeed (MachNumber (PSXC, PSXC*RRS), ATX)
# define a dataframe for the relative wind from the gust pod:
d <- data.frame("U_RW"=TAS_GP)
d["V_RW"] <- TAS_GP * tan (SS_GP * Cradeg)
d["W_RW"] <- TAS_GP * tan (AOA_GP * Cradeg)
d2 <- data.frame ("U2"=TASX)
d2["V2"] <- TASX * tan (SS * Cradeg)
d2["W2"] <- TASX * tan (AOA * Cradeg)



## ----Rotation-matrices---------------------------------------------------

rw <- as.matrix(d)
rw2 <- as.matrix(d2)
cosphi <- cos (CROLL_GP * Cradeg)
sinphi <- sin (CROLL_GP * Cradeg)
costheta <- cos (CPITCH_GP * Cradeg)
sintheta <- sin (CPITCH_GP * Cradeg)
cospsi <- cos (CTHDG_GP * Cradeg)
sinpsi <- sin (CTHDG_GP * Cradeg)
cosphi2 <- cos (ROLL * Cradeg)
sinphi2 <- sin (ROLL * Cradeg)
costheta2 <- cos (PITCH * Cradeg)
sintheta2 <- sin (PITCH * Cradeg)
cospsi2 <- cos (THDG * Cradeg)
sinpsi2 <- sin (THDG * Cradeg)
DL <- length(TASX)
One <- rep (1, DL)
Z <- rep (0, DL)
T1 <- array (c(One,Z,Z,Z,cosphi,-sinphi,Z,sinphi,cosphi), 
             dim=c(DL,3,3))
T2 <- array (c (costheta,Z,sintheta,Z,One,Z,-sintheta,Z,costheta),
             dim=c(DL,3,3))
T3 <- array (c (cospsi,-sinpsi,Z,sinpsi,cospsi,Z,Z,Z,One),
             dim=c(DL,3,3))
TT1 <- array (c(One,Z,Z,Z,cosphi2,-sinphi2,Z,sinphi2,cosphi2), 
             dim=c(DL,3,3))
TT2 <- array (c (costheta2,Z,sintheta2,Z,One,Z,-sintheta2,Z,costheta2),
             dim=c(DL,3,3))
TT3 <- array (c (cospsi2,-sinpsi2,Z,sinpsi2,cospsi2,Z,Z,Z,One),
             dim=c(DL,3,3))


## ----Calculate-new-wind-variables-WDG-WSG-WIG----------------------------

WDG <- vector ("numeric", DL)
WSG <- vector ("numeric", DL)
WIG <- vector ("numeric", DL)
WDX <- vector ("numeric", DL)
WSX <- vector ("numeric", DL)
WIX <- vector ("numeric", DL)
TASG <- vector ("numeric", DL)
Hlast <- 0.
# note: tried e.g. tensor (aperm(T1), rw, 2, 2) -- fails
#       because of size of vector that it tries to allocate
# This works, though:
# Y1t <- mapply(FUN="%*%", lapply(X=apply(aperm(T1), 3, as.data.frame), FUN=as.matrix, nrow=3, ncol=3), as.data.frame(aperm(rw)))
# etc.
for (i in 10000:10005) {
  Y1 <- aperm(T1[i,,]) %*% matrix (rw[i,], 3)
  Y2 <- aperm(T2[i,,]) %*% Y1
  Y3 <- aperm(T3[i,,]) %*% Y2
  WG <- matrix (c(-CVNS_GP[i], -CVEW_GP[i], CVSPD_GP[i]), 
                3, 1)
  Y4 <- Y3 + WG
  WDG[i] <- atan2 (Y4[2], Y4[1]) / Cradeg
  if ((!is.na(WDG[i])) & (WDG[i] < 0.)) {
    WDG[i] <- WDG[i] + 360.
  }
  WSG[i] <- sqrt (Y4[1]**2 + Y4[2]**2)
  WIG[i] <- Y4[3]
}



## ----Add-WIG-WDG-WSG-TASG-to-the-netCDF-file,echo=F,include=F------------
# add new variables to the netCDF file:

netCDFfile <- open.ncdf (fnew, write=TRUE)

if ("sps25" %in% names(netCDFfile$dim)) {
  dim(WIG) <- c(25, dim(TASX)/25)
  dim(WSG) <- c(25, dim(TASX)/25)
  dim(WDG) <- c(25, dim(TASX)/25)
  dim(TAS_GP) <- c(25, dim(TASX)/25)
  dim(WIX) <- c(25, dim(TASX)/25)
  dim(WSX) <- c(25, dim(TASX)/25)
  dim(WDX) <- c(25, dim(TASX)/25)

  varWIG <- var.def.ncdf ("WIG", "m/s", 
            c(netCDFfile$dim["sps25"], 
              netCDFfile$dim["Time"]), -32767.,
            "Vertical wind based on measurements from the gust pod")
  varWDG <- var.def.ncdf ("WDG", "deg.", 
            c(netCDFfile$dim["sps25"],
              netCDFfile$dim["Time"]), -32767.,
            "Wind direction based on measurements from the gust pod")
  varWSG <- var.def.ncdf ("WSG", "m/s", 
          c(netCDFfile$dim["sps25"],
            netCDFfile$dim["Time"]), -32767.,
          "Wind speed based on measurements from the gust pod")
  varWIX <- var.def.ncdf ("WIX", "m/s", 
            c(netCDFfile$dim["sps25"], 
              netCDFfile$dim["Time"]), -32767.,
            "Vertical wind based on measurements from the radome, revised")
  varWDX <- var.def.ncdf ("WDX", "deg.", 
            c(netCDFfile$dim["sps25"],
              netCDFfile$dim["Time"]), -32767.,
            "Wind direction based on measurements from the radome, revised")
  varWSX <- var.def.ncdf ("WSX", "m/s", 
          c(netCDFfile$dim["sps25"],
            netCDFfile$dim["Time"]), -32767.,
          "Wind speed based on measurements from the radome, revised")
  varTASG <- var.def.ncdf ("TASG", "m/s", 
          c(netCDFfile$dim["sps25"],
            netCDFfile$dim["Time"]), -32767.,
          "True Airspeed, gust pod")
} else {
  varWIG <- var.def.ncdf ("WIG", "m/s", 
            netCDFfile$dim["Time"], -32767.,
            "Vertical wind based on measurements from the gust pod")
  varWDG <- var.def.ncdf ("WDG", "deg.", 
            netCDFfile$dim["Time"], -32767.,
            "Wind direction based on measurements from the gust pod")
  varWSG <- var.def.ncdf ("WSG", "m/s", 
          netCDFfile$dim["Time"], -32767.,
          "Wind speed based on measurements from the gust pod") 
  varWIX <- var.def.ncdf ("WIX", "m/s", 
            netCDFfile$dim["Time"], -32767.,
            "Vertical wind based on measurements from the radome, revised")
  varWDX <- var.def.ncdf ("WDX", "deg.", 
            netCDFfile$dim["Time"], -32767.,
            "Wind direction based on measurements from the radome, revised")
  varWSX <- var.def.ncdf ("WSX", "m/s", 
          netCDFfile$dim["Time"], -32767.,
          "Wind speed based on measurements from the radome, revised") 
  varTASG <- var.def.ncdf ("TASG", "m/s", 
          netCDFfile$dim["Time"], -32767.,
          "True Airspeed, gust pod") 
}
newfile <- var.add.ncdf (netCDFfile, varWIG)
newfile <- var.add.ncdf (newfile, varWDG)
newfile <- var.add.ncdf (newfile, varWSG)
newfile <- var.add.ncdf (newfile, varWIX)
newfile <- var.add.ncdf (newfile, varWDX)
newfile <- var.add.ncdf (newfile, varWSX)
newfile <- var.add.ncdf (newfile, varTASG)
put.var.ncdf (newfile, "WIG", WIG)
put.var.ncdf (newfile, "WDG", WDG)
put.var.ncdf (newfile, "WSG", WSG)
put.var.ncdf (newfile, "WIX", WIX)
put.var.ncdf (newfile, "WDX", WDX)
put.var.ncdf (newfile, "WSX", WSX)
put.var.ncdf (newfile, "TASG", TAS_GP)
close.ncdf (newfile)

detach(Data)
rm(WIG, WDG, WSG, WIX, WDX, WSX, TAS_GP)



## ----Plotted-Results, fig.cap=c("WIG (vertical wind based on the gust pod) plotted against WIC (vertical wind from the conventional radome-based gust system)", "Angle of attack determined from gust-pod measurements, plotted vs. corresponding measurements AKRD from the standard wind sensing system", "Comparison of horizontal wind measurements from the gust pod (red lines) and from the standard wind measurements WDC and WSC (thicker blue lines).", "Comparison of vertical wind measurements from the gust pod (red line) and from the standard wind measurement (WIC, blue line)."), out.width="300pt", fig.align="center",echo=FALSE----

print (c("Flight processed in this run: ",Flight))
Data <- getNetCDF (fnew, c(VarNames,"WIG", "WDG", "WSG", "TASG", "WIX", "WDX", "WSX"))
Valid <- (Data$TASX > 130.)
plot (Data$WIC[Valid], Data$WIG[Valid], ylab="WIG (gust pod vertical wind, m/s)",
      xlim=c(-5.,5.), ylim=c(-5.,5.),
      xlab="WIC, standard vertical wind [m/s]", pch=16, 
      cex=0.5, col='blue')
lines (c(-5.,5.), c(-5.,5.), lty=2, lwd=2, col='orange')
print (mean (Data$WIG[Valid]-Data$WIC[Valid], na.rm=TRUE))
print (sd (Data$WIG[Valid]-Data$WIC[Valid], na.rm=TRUE))
# text (1.,-3.5, "mean WIG-WIC: -0.21 m/s\nst dev 0.25 m/s")
plot (Data$AKRD[Valid], AOA_GP[Valid], pch=16, cex=0.5, xlab="AKRD [deg.]", 
      xlim=c(1.,5.), ylim=c(-3.,2.),
      ylab="AOA from GP [deg.]")
lines (c(2.0, 3.6), c(2.0-3.65, 3.6-3.66), lty=2, lwd=2, 
       col='orange')

op <- par (mfrow=c(2,1), mar=c(4,4,0,2)+0.1)
plotWAC (Data$Time, Data$WDC, ylab="Wind Direction [deg.]")
points (Data$Time, Data$WDG, type='l', col='red')
legend ("bottomright", legend=c("WDC", "WDG"), lty=1, col=c('blue', 'red'))
plotWAC (Data$Time, Data$WSC, ylab="Wind Speed [m/s]")
points (Data$Time, Data$WSG, type='l', col='red')
legend ("bottomright", legend=c("WSC", "WSG"), lty=1, col=c('blue', 'red'))
op <- par (mfrow=c(1,1), mar=c(4,4,5,2)+0.1)
plotWAC (Data$Time, Data$WIC, ylab="Vertical Wind [m/s]")
points (Data$Time, Data$WIG, type='l', col='red')
legend ("bottomright", legend=c("WIC", "WIG"), lty=1, col=c('blue', 'red'))
# sd (Data$WSC-WSG, na.rm=TRUE)
# sd (Data$WDC-WDG, na.rm=TRUE)



