require(Ranadu)

Flight <- "rf21"
SavePlotsToFiles <- FALSE

fname = sprintf("/Data/DEEPWAVE/DEEPWAVE%s.nc", Flight)
VarNames <- standardVariables()
VarNames <- c(VarNames, "ATHR1", "ATHR2", "ATRL", "AT_A", "AT_A2")
VarNames <- c(VarNames, "EW_DPL", "DP_DPL", "PALTF", "EW_DPR",
              "DP_DPR", "EW_VXL", "DP_VXL", "CAVP_CR2", 
              "MIRRTMP_CR2", "DP_CR2")
VarNames <- c(VarNames, "PSFC", "PS_A", "PS_GP", "QCFC", 
              "QCRC", "QC_A", "QC_GP")
VarNames <- c(VarNames, "GGQUAL", "GGVEW", "GGVNS", "VEW", 
              "VNS", "AKRD", "SSLIP")
VarNames <- c(VarNames, "CONCU_RWO", "CONCU100_RWO", "CONCU100_RWO", 
              "CONCU500_RWO", "CNTS")
VarNames <- c(VarNames, "CLAT_LAMS", "CLON_LAMS", "CPITCH_LAMS",
              "CROLL_LAMS", "CTHDG_LAMS", "CVEW_LAMS", "CVNS_LAMS",
              "CVSPD_LAMS", "PITCH", "ROLL", "THDG")
Data <- getNetCDF (fname, VarNames)
attach(Data)
#netCDFfile = open.ncdf(fname)
#V <- get.var.ncdf(netCDFfile, "BEAM1_LAMS")
# construct flight track with map
require(maps)
require(mapdata)
require(mapproj)
if (SavePlotsToFiles) {
  png (filename=sprintf("~/RStudio/DEEPWAVE/DW%sPlot%%02d.png", Flight))
}
op <- par (mfrow=c(1,1), mar=c(5,5,2,2)+0.1)
plotTrack (LONC, LATC, Time, Spacing=60, WindFlags=2)
title (Flight)
#plotTrack (LONC, LATC, Time, setRange (Time, 25000, 43000), 
#           xc=-112., yc=42.0, sz=4., WindFlags=3)
# check temperatures
plotWAC(Time, ATHR1, col="blue")
points (Time, ATHR2, type='l', col="black")
points (Time, ATRL, type='l', col="red")
points (Time, AT_A, type='l', col="green")
legend (median(Time), max(ATHR1, na.rm=TRUE)-10., 
        legend=c("ATHR1", "ATHR2", "ATRL", "AT_A"), lty=1,
        col=c("blue", "black", "red", "green"))
# eliminate points with TASX<100.
ATHR1 <- ATHR1[TASX > 100.]
ATHR2 <- ATHR2[TASX > 100.]
ATRL <- ATRL[TASX > 100.]
AT_A <- AT_A[TASX > 100.]
#ATHR1,2 section
plot (ATHR1, ATHR2, pch=16, cex=0.5)
lines (c(-70.,30.), c(-70.,30.), col="orange")
points (ATHR1, (ATHR2-ATHR1)*10, pch=16, cex=0.5, col='red')
lines (c(-70.,30.), c(5.,5.), lty=2, col='red')
lines (c(-70.,30.), c(-5.,-5.), lty=2, col='red')
legend ("bottomright",
        legend=c("(ATHR2-ATHR1)*10", "+/-0.5"), lty=c(1,2),
        col=c("red", "red"), cex=0.7)
fm <- lm(ATHR2~ATHR1)
coef <- coefficients (fm)
if (coef[1] < 0.) {
  t <- sprintf ("ATHR2=%.3f(ATHR1)%.3f\nmean diff %.2f +/- %.2f", coef[2], 
                coef[1], mean (ATHR2-ATHR1, na.rm=TRUE),
                sd(ATHR2-ATHR1, na.rm=TRUE))
} else {
  t <- sprintf ("ATHR2=%.3f(ATHR1)+%.3f\nmean diff %.2f +/-%.2f", coef[2], 
                coef[1], mean (ATHR2-ATHR1, na.rm=TRUE),
                sd(ATHR2-ATHR1, na.rm=TRUE))
}
title(t)

#ATRL section:
plot (ATHR1, ATRL, pch=16, cex=0.5)
lines (c(-70.,30.), c(-70.,30.), col="orange")
points (ATHR1, (ATRL-ATHR1)*10, pch=16, cex=0.5, col='red')
lines (c(-70.,30.), c(5.,5.), lty=2, col='red')
lines (c(-70.,30.), c(-5.,-5.), lty=2, col='red')
legend ("bottomright", cex=0.6,
        legend=c("(ATRL-ATHR1)*10", "+/-0.5"), lty=c(1,2),
        col=c("red", "red"))
fm <- lm(ATRL~ATHR1)
coef <- coefficients (fm)
if (coef[1] < 0.) {
  t <- sprintf ("ATRL=%.3f(ATHR1)%.3f\n mean diff %.2f +/-%.2f", coef[2], 
                coef[1], mean(ATRL-ATHR1, na.rm=TRUE),
                sd(ATRL-ATHR1, na.rm=TRUE))
} else {
  t <- sprintf ("ATRL=%.3f(ATHR1)+%.3f\n mean diff %.2f +/-%.2f", coef[2], 
                coef[1], mean(ATRL-ATHR1, na.rm=TRUE), 
                sd(ATRL-ATHR1, na.rm=TRUE))
}
title(t)

#AT_A section:
plot (ATHR1, AT_A, pch=16, cex=0.5)
lines (c(-70.,30.), c(-70.,30.), col="orange")
points (ATHR1, (AT_A-ATHR1)*10, pch=16, cex=0.5, col='red')
lines (c(-70.,30.), c(5.,5.), lty=2, col='red')
lines (c(-70.,30.), c(-5.,-5.), lty=2, col='red')
legend ("bottomright", 
        legend=c("(AT_A-ATHR1)*10", "+/-0.5"), lty=c(1,2),
        cex=0.6, col=c("red", "red"))
fm <- lm(AT_A~ATHR1)
coef <- coefficients (fm)
if (coef[1] < 0.) {
  t <- sprintf ("AT_A=%.3f(ATHR1)%.3f\n mean diff %.2f +/-%.2f", coef[2], 
                coef[1], mean(AT_A-ATHR1, na.rm=TRUE),
                sd(AT_A-ATHR1, na.rm=TRUE))
} else {
  t <- sprintf ("AT_A=%.3f(ATHR1)+%.3f\n mean diff %.2f +/-%.2f", coef[2], 
                coef[1], mean(AT_A-ATHR1, na.rm=TRUE), 
                sd(AT_A-ATHR1, na.rm=TRUE))
}
title(t)

#humidity
plotWAC(Time, DP_VXL)
points (Time, DP_DPL, type='l', col='green')
points (Time, DP_DPR, type='l', col='cyan')
points (Time, DP_CR2, type='l', col='red')
legend ("top",
        legend=c("DP_VXL", "DP_DPL", "DP_DPR", "DP_CR2"),
        lty=1, cex=0.6,
        col=c("blue", "green", "cyan", "red"))

#pressure
plotWAC (Time, PSFC)
points (Time, PS_A, type='l', col='green')
points (Time, PS_GP, type='l', col='cyan')
points (Time, (PSFC-PS_A)*100+600, type='l', col='red')
lines (c(Time[1],Time[length(Time)]), c(600.+100., 600.+100.), type='l', lty=2, col='red')
lines (c(Time[1],Time[length(Time)]), c(600.-100., 600.-100.), type='l', lty=2, col='red')
legend ("top", legend=c("PSFC", "PS_A", "PS_GP", "(PSFC-PS_A)*100+600", "+/-1 hPa"),
        lty=c(1,1,1,1,2), cex=0.55,
        col=c('blue', 'green', 'cyan', 'red', 'red'))

#dynamic pressure
plotWAC (Time, QCFC)
points (Time, QC_A, type='l', col='green')
points (Time, QC_GP, type='l', col='cyan')
points (Time, QCRC, type='l', col='orange')
points (Time, (QCFC-QC_A)*10+60, type='l', col='red')
lines (c(Time[1],Time[length(Time)]), c(60.+10., 60.+10.), type='l', lty=2, col='red')
lines (c(Time[1],Time[length(Time)]), c(60.-10., 60.-10.), type='l', lty=2, col='red')
legend ("top", legend=c("QCFC", "QC_A", "QC_GP",
                        "QCRC", "(QCFC-QC_A)*10+60", "+/-1 hPa"),
        lty=c(1,1,1,1,1,2), cex=0.55,
        col=c('blue', 'green', 'cyan', 'orange', 'red', 'red'))

#wind
op <- par (mfrow=c(3,1), mar=c(4,4,0,2)+0.1)
plotWAC (Time, WDC, ylab="WDC [deg.]")
plotWAC (Time, WSC, ylab="WSC [m/s]")
plotWAC (Time, WIC, ylab="WIC [m/s]")
op <- par (mfrow=c(2,1), mar=c(4,4,0,2)+0.1)
plotWAC (Time, GGVEW, ylab="VEW")
points (Time, VEW, type='l', col='green', lty=2, lwd=2)
points (Time, (GGVEW-VEW)*50, type='l', col='red')
lines (c(Time[1], Time[length(Time)]), c(50.,50.), type='l', lty=2)
lines (c(Time[1], Time[length(Time)]), c(-50.,-50.), type='l', lty=2)
plotWAC (Time, GGVNS, ylab="VNS")
points (Time, VNS, type='l', col='green', lty=2, lwd=2)
points (Time, (GGVNS-VNS)*50, type='l', col='red')
lines (c(Time[1], Time[length(Time)]), c(50.,50.), type='l', lty=2)
lines (c(Time[1], Time[length(Time)]), c(-50.,-50.), type='l', lty=2)
op <- par (mfrow=c(3,1), mar=c(4,4,0,2)+0.1)
plotWAC (Time, GGQUAL, ylab='GGQUAL')
plotWAC (Time, AKRD, ylab="AKRD [deg.]")
plotWAC (Time, SSLIP, ylab="SSLIP [deg.]")

op <- par (mfrow=c(3,1), mar=c(4,4,0,2)+0.1)
plot (Time, PITCH, type='l', ylim=c(-10.,10.), col="blue", ylab="PITCH [deg.]")
points (Time, CPITCH_LAMS, type='l', lty=2, col="green")
#points (Time, THDG/18.-10., type='l', col="red")
plotWAC (Time, ROLL, yl="ROLL [deg.]")
points (Time, CROLL_LAMS, type='l', lty=2, col="green" )
plotWAC (Time, THDG, yl="Heading")
points (Time, CTHDG_LAMS, type='l', lty=2, col="green")
legend ("topright", legend=c("THDG", "CTHDG"), lty=c(1,2), col=c("blue", "green"))
#legend (Time[getIndex(Time, 120000)], 450., legend=c("THDG", "CTHDG"), lty=c(1,2), col=c("blue", "green"), xpd=NA, horiz=TRUE)

#UHSAS
op <- par (mfrow=c(1,1), mar=c(4,5,2,2)+0.1)
plot (Time, CONCU_RWO, ylab="CONCU", type='l', log='y')
points (Time, CONCU100_RWO, type='l', col="green")
points (Time, CONCU500_RWO, type='l', col="red")

#CN
plot (Time, CNTS, ylab="CN Concentration", type='l', log='y')
# netCDFfile = open.ncdf(fname)
# CUHSAS <- get.var.ncdf(netCDFfile, "CUHSAS_RWO")
# CellSizes <- att.get.ncdf (netCDFfile, "CUHSAS_RWO", "CellSizes")
# CellLimits <- CellSizes$value
# 
# 
# for (j in seq(1400,1420,2)) {
#   UHSAS <- vector ("numeric", 100)
#   for (k in 0:60) {
#     UHSAS <- UHSAS + CUHSAS[,getIndex(Time, j*100+k)]
#   }
#   UHSAS <- UHSAS / 60.
#   UHSAS[UHSAS <= 0.] <- 1.e-3
#   plot (CellLimits, UHSAS, type='l', ylim=c(1.e-3,1.e1), xlab="Diameter [um]", log="xy")
#   lines (c(0.5,0.5), c(1.e-4, 1.e3), lty=2, col="darkgreen")
#   lines (c(0.1,0.1), c(1.e-4, 1.e3), lty=2, col="darkgreen")
#   title (sprintf ("1-min average starting at %d", j))
# }
# close.ncdf (netCDFfile)

if (SavePlotsToFiles) {
  dev.off()
}
detach (Data)
rm (AT_A, ATHR1, ATHR2, ATRL)



