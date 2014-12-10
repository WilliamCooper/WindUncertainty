# recover angle of attack from gust-probe measurements?

# AlgorithmNotes ----------------------------------------------------------


require(Ranadu)
Flight <- "rf05"
SavePlotsToFiles <- FALSE
fname = sprintf("/home/Data/DEEPWAVE/DEEPWAVE%s.nc", Flight)
VarNames <- standardVariables()
VarNames <- c(VarNames, "ATHR1", "ATHR2", "ATRL", "AT_A", "AT_A2")
VarNames <- c(VarNames, "EW_DPL", "DP_DPL", "PALTF", "EW_DPR",
              "DP_DPR", "EW_VXL", "DP_VXL", "CAVP_CR2", 
              "MIRRTMP_CR2")
VarNames <- c(VarNames, "PSFC", "PS_A", "PS_GP", "QCFC", 
              "QCRC", "QC_A", "QC_GP")
VarNames <- c(VarNames, "GGQUAL", "GGVEW", "GGVNS", "VEW", 
              "VNS", "AKRD", "SSLIP")
VarNames <- c(VarNames, "CONCU_RWO", "CONCU100_RWO", "CONCU100_RWO", 
              "CONCU500_RWO", "CNTS")
VarNames <- c(VarNames, "CLAT_LAMS", "CLON_LAMS", "CPITCH_LAMS",
              "CROLL_LAMS", "CTHDG_LAMS", "CVEW_LAMS", "CVNS_LAMS",
              "CVSPD_LAMS", "PITCH", "ROLL", "THDG")
VarNames <- c(VarNames, "ADIFR_GP", "BDIFR_GP", "PS_GP", 
              "CROLL_GP", "CTHDG_GP", "CVNS_GP", "CVEW_GP",
              "QC_GP", "WI_GP", "CPITCH_GP", "VSPD", "CVSPD_GP")
Data <- getNetCDF (fname, VarNames)
i1 <- getIndex(Data$Time, 74600)
# Valid <- (Data$Time > Data$Time[i1])
# Valid <- Valid & ((Data$TASX > 130.) & (abs (Data$ROLL) < 10.) &
#             (abs(Data$VSPD) < 10))
Valid <- ((Data$TASX > 130.) & (abs (Data$ROLL) < 10.) &
                 (abs(Data$VSPD) < 2))
D <- Data[Valid,]
attach(D)
Cradeg <- pi / 180.
AQR <- ADIFR_GP/QC_GP
RR2 <- QC_GP/PS_GP
Mach <- MachNumber (PSXC, QCXC)
fm <- lm (AKRD~AQR+RR2+PS_GP+Mach+QC_GP)
cf2 <- coefficients (fm)
A <- cf2[1]+cf2[2]*AQR+cf2[3]*RR2+cf2[4]*PS_GP+cf2[5]*Mach+cf2[6]*QC_GP
plot (AKRD, A, pch=16, cex=0.5, col='blue')
WIX <- TASX*(A-PITCH)*Cradeg + VSPD
plot (WIC, WI_GP, pch=16, col='green', cex=0.5)
lines (c(-5.,5.), c(-5.,5.), lty=2, col='orange', lw=2)
points (WIC, WIX, col='blue', pch=16, cex=0.5)
legend("bottomright", legend=c("WIX -- new calc", "WI_GP", "1:1 reference line"),
       pch=c(16,16,1), pt.cex=c(0.5,0.5,0.0), lty=c(0, 0, 2), 
       col=c("blue", "green", "orange"), lwd=c(0,0,2))

AOAREF <- CPITCH_GP - asin(CVSPD_GP/TASX) / Cradeg
plot (AKRD, AOAREF, pch=16, cex=0.5, col='blue')

fm <- lm (AOAREF~AQR+RR2+PS_GP+QC_GP+Mach)
cf3 <- coefficients (fm)
A2 <- cf3[1]+cf3[2]*AQR+cf3[3]*RR2+cf3[4]*PS_GP+cf3[5]*QC_GP+cf3[6]*Mach
plot (AOAREF, A2, pch=16, cex=0.5, col='blue')
lines (c(-3.,3.), c(-3.,3.), lty=2, lwd=2, col='orange')
fm3 <- lm(AOAREF~AQR)
cf <- coefficients (fm3)
A3 <- cf[1]+cf[2]*AQR
plot (AOAREF, A3, pch=16, cex=0.5, col='blue')
lines (c(-3.,3.), c(-3.,3.), lty=2, lwd=2, col='orange')
title ("DEEPWAVE rf05\nStandard error 0.13 deg")
legend ("bottomright", legend="1:1 reference line", lty=2, col='orange', lwd=2)
B1 <- QCXC/PSXC
B2 <- QC_GP/PS_GP
fm4 <- lm(B1~B2+I(B2^2)+AQR)
fm4 <- lm(B1~B2+I(B2^2)+AQR)
cf4 <- coefficients (fm4)
B3 <- cf4[1]+cf4[2]*B2+cf4[3]*B2**2+cf4[4]*AQR
plot (B1, B3, pch=16, cex=0.5, col='blue')
lines (c(0., 1.), c(0.,1.), lwd=2, lty=2, col='orange')
Vg <- TrueAirspeed (MachNumber (PSXC, PSXC*B3), ATX)
plot(TASX, Vg, pch=16, cex=0.5, col='blue')
lines (c(100.,300.), c(100.,300.), lty=2, col='orange')
sd(TASX-Vg)
sd((TASX-Vg)[abs(VSPD) < 2])

# now get the relative wind components
detach(D)
attach(Data)
AQR <- ADIFR_GP/QC_GP
RR2 <- QC_GP/PS_GP
SS <- cf[2] * BDIFR_GP / QC_GP-2.9
B2 <- QC_GP/PS_GP
B3 <- cf4[1]+cf4[2]*B2+cf4[3]*B2**2+cf4[4]*AQR
Vg <- TrueAirspeed (MachNumber (PSXC, PSXC*B3), ATX)
Mach <- MachNumber (PSXC, QCXC)
A2 <- cf3[1]+cf3[2]*AQR+cf3[3]*RR2+cf3[4]*PS_GP+cf3[5]*QC_GP+cf3[6]*Mach

A3 <- cf[1]+cf[2]*AQR
d <- data.frame("U_RW"=Vg)
d["V_RW"] <- Vg * tan(SS * Cradeg)
d["W_RW"] <- Vg * tan (A3 * Cradeg)
rw <- as.matrix(d)
cosphi <- cos (CROLL_GP * Cradeg)
sinphi <- sin (CROLL_GP * Cradeg)
costheta <- cos (CPITCH_GP * Cradeg)
sintheta <- sin (CPITCH_GP * Cradeg)
cospsi <- cos (CTHDG_GP * Cradeg)
sinpsi <- sin (CTHDG_GP * Cradeg)
DL <- length(TASX)
One <- rep (1, DL)
Z <- rep (0, DL)
T1 <- array (c(One,Z,Z,Z,cosphi,-sinphi,Z,sinphi,cosphi), 
             dim=c(DL,3,3))
T2 <- array (c (costheta,Z,sintheta,Z,One,Z,-sintheta,Z,costheta),
             dim=c(DL,3,3))
T3 <- array (c (cospsi,-sinpsi,Z,sinpsi,cospsi,Z,Z,Z,One),
             dim=c(DL,3,3))
i <- 8000
WDG <- vector ("numeric", DL)
WSG <- vector ("numeric", DL)
WIG <- vector ("numeric", DL)
Hlast <- 0.
for (i in 1:DL) {
  rwx <- matrix(rw[i,], 3, 1)
  Y1 <- aperm(T1[i,,]) %*% matrix (rw[i,], 3, 1)
  Y2 <- aperm(T2[i,,]) %*% Y1
  Y3 <- aperm(T3[i,,]) %*% Y2
#Y3 <- T3[i,,] %*% Y2
#print (Y3)
  WG <- matrix (c(-CVNS_GP[i], -CVEW_GP[i], CVSPD_GP[i]), 
                3, 1)
  Y4 <- Y3 + WG
#print (Y4)
  WDG[i] <- atan2 (Y4[2], Y4[1]) / Cradeg
  if ((!is.na(WDG[i])) & (WDG[i] < 0.)) {
    WDG[i] <- WDG[i] + 360.
  }
  WSG[i] <- sqrt (Y4[1]**2 + Y4[2]**2)
  WIG[i] <- Y4[3]
  if (is.na(CTHDG_GP[i]) | ((abs(CTHDG_GP[i]-Hlast) > 3.) 
                             & (abs(Hlast-180.) < 5.))) {
    WIG[i] <- NA
    WDG[i] <- NA
    WSG[i] <- NA
  }
  Hlast <- CTHDG_GP[i]
#print (c(WD, WS))
#print (c(D$WDC[i], D$WSC[i]))
}

plot (Data$WIC, WIG, pch=16, cex=0.5, col='blue')
lines (c(-5.,5.), c(-5.,5.), lty=2, lwd=2, col='orange')
print (mean (WIG-Data$WIC, na.rm=TRUE))
print (sd (WIG-Data$WIC, na.rm=TRUE))
text (1.,-3.5, "mean WIG-WIC: -0.21 m/s\nst dev 0.25 m/s")
plot (Data$AKRD, A2, pch=16, cex=0.5, xlab="AKRD [deg.]", 
      ylab="AOA from GP [deg.]")
lines (c(2.0, 3.6), c(2.0-3.65, 3.6-3.66), lty=2, lwd=2, 
       col='orange')
op <- par (mfrow=c(2,1), mar=c(4,4,0,2)+0.1)
plotWAC (Data$Time, Data$WDC, ylab="Wind Direction [deg.]")
points (Data$Time, WDG, type='l', col='red')
legend ("bottomright", legend=c("WDC", "WDG"), lty=1, col=c('blue', 'red'))
plotWAC (Data$Time, Data$WSC, ylab="Wind Speed [m/s]")
points (Data$Time, WSG, type='l', col='red')
legend ("bottomright", legend=c("WSC", "WSG"), lty=1, col=c('blue', 'red'))
op <- par (mfrow=c(2,1), mar=c(4,4,0,2)+0.1)
plotWAC (Data$Time, Data$WDC, ylab="Wind Direction [deg.]")
points (Data$Time, WDG, type='l', col='red')
legend ("bottomright", legend=c("WDC", "WDG"), lty=1, col=c('blue', 'red'))
plotWAC (Data$Time, Data$WSC, ylab="Wind Speed [m/s]")
points (Data$Time, WSG, type='l', col='red')
legend ("bottomright", legend=c("WSC", "WSG"), lty=1, col=c('blue', 'red'))
sd (Data$WSC-WSG, na.rm=TRUE)
sd (Data$WDC-WDG, na.rm=TRUE)

# add new variables to the netCDF file:
AQR <- ADIFR_GP/QC_GP
RR2 <- QC_GP/PS_GP
Mach <- MachNumber (PSXC, QCXC)
# use coefficients from above
A <- cf[1]+cf[2]*AQR+cf[3]*RR2+cf[4]*PS_GP+cf[5]*Mach+cf[6]*QC_GP
WIX <- TASX*(A-PITCH)*Cradeg + VSPD
netCDFfile <- open.ncdf ("/home/Data/DEEPWAVE/NEWrf05.nc", write=TRUE)
#Time <- get.var.ncdf (netCDFfile, "Time")
varWIG <- var.def.ncdf ("WIG", "m/s", netCDFfile$dim["Time"], -32767., "Vertical wind based on measurements from the gust pod")
varWDG <- var.def.ncdf ("WDG", "deg.", netCDFfile$dim["Time"], -32767., "Wind direction based on measurements from the gust pod")
varWSG <- var.def.ncdf ("WSG", "m/s", netCDFfile$dim["Time"], -32767., "Wind speed based on measurements from the gust pod")
newfile <- var.add.ncdf (netCDFfile, varWIG)
newfile <- var.add.ncdf (newfile, varWDG)
newfile <- var.add.ncdf (newfile, varWSG)
#close.ncdf (newfile)
#newfile <- open.ncdf ("/home/Data/DEEPWAVE/NEWrf05.nc", write=TRUE)
put.var.ncdf (newfile, "WIG", WIG, count=-1)
put.var.ncdf (newfile, "WDG", WDG)
put.var.ncdf (newfile, "WSG", WSG)
close.ncdf (newfile)

#ncNew <- create.ncdf ("/home/Data/DEEPWAVE/NEWGrf05.nc", list(varWIG, varWDG,varWSG,netCDFfile))
#plot (d$U_RW)


detach(Data)

