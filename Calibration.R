require(Ranadu)

Flight <- "rf15"
SavePlotsToFiles <- TRUE
fname = sprintf("/home/Data/DEEPWAVE/DEEPWAVE%s.nc", Flight)
VarNames <- standardVariables()
VarNames <- c("VYC","GGALT","LATC", "LONC", "PSXC", "QCXC",
              "WDC", "WSC", "GGVEW", "GGVNS", "VEW", "VNS", "TASX",
              "ADIFR", "AKRD", "SSLIP", "PITCH", 
              "ROLL", "THDG", "BDIFR",
              "ADIFR_GP", "BDIFR_GP", "PS_GP", "QC_GP",
              "CROLL_GP", "CPITCH_GP", "CTHDG_GP", "WIC",
              "CVNS_GP", "CVEW_GP", "VSPD", "CVSPD_GP",
              "ATX")
D <- getNetCDF (fname, VarNames)
r1 <- setRange (D$Time, 32100,32900)
r2 <- setRange (D$Time, 41500,42300)
r3 <- setRange (D$Time, 50100,51100)
r <- c(r1,r2,r3)
Data <- D[r,]
Valid <- (Data$TASX > 130.)
Data <- Data[Valid,]
attach(Data)
Cradeg <- pi / 180.
Mach <- MachNumber (PSXC, QCXC)  # uses conventional q, p
AOAREF <- PITCH - asin(VSPD/TASX) / Cradeg
AQR <- ADIFR/QCXC # basic pressure ratio for AOA
RR2 <- QCXC/PSXC    # q/p ratio
plot (AOAREF, AKRD, pch=16, cex=0.5, col='blue')
#fm <- lm (AOAREF~AQR+RR2+PS_GP+QC_GP+Mach)
fm <- lm (AOAREF~AQR+RR2+Mach)
cf3 <- coefficients (fm)
A2 <- cf3[1]+cf3[2]*AQR+cf3[3]*RR2+cf3[4]*Mach
plotWAC (AOAREF, A2, xlab="Predicted AOA [deg.]", 
         ylab="Fit AOA [deg.]", type='p', pch=16, cex=0.5, col='blue')
lines (c(-3.,6.), c(-3.,6.), lty=2, lwd=2, col='orange')
title("Radome")
AOAREF_GP <- CPITCH_GP - asin(CVSPD_GP/TASX) / Cradeg
AQR_GP <- ADIFR_GP/QC_GP # basic pressure ratio for AOA
RR2_GP <- QC_GP/PS_GP    # q/p ratio from the gust pod
#fm <- lm (AOAREF~AQR+RR2+PS_GP+QC_GP+Mach)
fm_GP <- lm (AOAREF_GP~AQR_GP+RR2_GP+Mach)
cf3_GP <- coefficients (fm_GP)
A2_GP <- cf3_GP[1]+cf3_GP[2]*AQR_GP+cf3_GP[3]*RR2_GP+cf3_GP[4]*Mach
plotWAC (AOAREF_GP, A2_GP, xlab="Predicted AOA [deg.]", 
         ylab="Fit AOA [deg.]", type='p', pch=16, cex=0.5, col='blue')
lines (c(-3.,6.), c(-3.,6.), lty=2, lwd=2, col='orange')
title ("Gust Pod")

# sideslip calibration:

detach(Data)
r1 <- setRange (D$Time, 33230,33500)
r2 <- setRange (D$Time, 43100,43300)
r3 <- setRange (D$Time, 52640,52930)
r <- c(r1,r2,r3)
#r <- r1
Data <- D[r,]
Valid <- (Data$TASX > 130.)
Data <- Data[Valid,]
attach(Data)
Hmean <- mean(THDG)
u <- -1. * WSC * sin (WDC*Cradeg)
v <- -1. * WSC * cos (WDC*Cradeg)
SSREF <- -THDG + atan2((VEW-u), (VNS-v))/ Cradeg
SSREF[SSREF < -180.] <- SSREF[SSREF < -180.] + 360.
plot (SSREF, SSLIP, pch=16, cex=0.6,col='blue')
BQR <- BDIFR / QCXC
sfm <- lm(SSREF~BQR)
summary(sfm)
cfs <- coefficients(sfm)
S2 <- cfs[1] + cfs[2] * BQR
plot (SSREF, S2, pch=16, cex=0.6, col='blue',
      xlab="SS Reference", ylab="SSLIP from fit")
lines(c(-2.,2.), c(-2.,2.), col='orange', lty=2, lwd=2)
title ("Radome")

SSREF_GP <- -CTHDG_GP + atan2((CVEW_GP-u), (CVNS_GP-v))/ Cradeg
SSREF_GP[SSREF_GP < -180.] <- SSREF_GP[SSREF_GP < -180.] + 360.
plot (SSREF_GP, SSLIP, pch=16, cex=0.6,col='blue')
BQR_GP <- BDIFR_GP / QC_GP
Mach <- MachNumber (PSXC, QCXC)
#sfm_GP <- lm(SSREF_GP~BQR_GP)
RR2_GP <- QC_GP / PS_GP
sfm_GP <- lm(SSREF_GP~BQR_GP)
summary(sfm_GP)
cfs_GP <- coefficients(sfm_GP)
S2_GP <- cfs_GP[1] + cfs_GP[2] * BQR_GP #+ cfs_GP[3] * RR2_GP + cfs_GP[4] * Mach
plot (SSREF_GP, S2_GP, pch=16, cex=0.6, col='blue',
      xlab="SS Reference", ylab="SSLIP from fit")
lines(c(-2.,2.), c(-2.,2.), col='orange', lty=2, lwd=2)
title ("Gust Pod")
detach(Data)
