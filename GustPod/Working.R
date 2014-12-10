require(Ranadu) 
library(knitr) 
require(ggplot2) 
require(grid) 
require(ggthemes) 
require(vioplot) 
require(plyr) 
require(signal) 
require(stats)
# set global chunk options 
#opts_chunk$set(fig.path='figure/WU-', fig.align='center', fig.show='hold')
#options(replace.assign=TRUE,width=49, digits=4)
Flight <- "rf15"        # this was the flight with cal maneuvers
Directory = "/home/Data/"
fname = sprintf("%sDEEPWAVE/DW%s_GPx.nc", Directory, Flight)
VarNames <- c("VYC","GGALT","LATC", "LONC", "PSXC", "QCXC",
              "WDC", "WSC", "GGVEW", "GGVNS", "VEW", "VNS", "TASX",
              "ADIFR", "AKRD", "SSLIP", "PITCH", "PSF", "QCF",
              "ROLL", "THDG", "BDIFR", "EWX", "GGVSPDB",
              "ADIFR_GP", "BDIFR_GP", "PS_GP", "QC_GP",
              "CROLL_GP", "CPITCH_GP", "CTHDG_GP", "WIC",
              "CVNS_GP", "CVEW_GP", "VSPD", "CVSPD_GP",
              "ATX", "GGVEWB", "GGVNSB", "VNSC", "VEWC",
              "WDX", "WSX", "WIX")
D <- getNetCDF (fname, VarNames, F=15)      # this handles high-rate data also
r <- setRange(D$Time, 33800,35500)
D <- D[r,]
cfs <- c(0.016, 21.9915)
SS <- cfs[1] + cfs[2]*BDIFR/QCF
plotWAC(D$Time, SS, ylab="Sideslip")
lines(c(D$Time[1], D$Time[length(D$Time)]), c(-0.02,-0.02), col='darkorange', lty=2, lwd=2)
yr <- mean(SS[setRange(D$Time,34000,34500)])
lines(c(D$Time[getIndex(D$Time,34000)], D$Time[getIndex(D$Time, 34500)]), c(yr,yr), col='darkorange', lty=2, lwd=4)
yr <- mean(SS[setRange(D$Time, 34730,35330)])
lines(c(D$Time[getIndex(D$Time,34730)], D$Time[getIndex(D$Time, 35330)]), c(yr,yr), col='darkorange', lty=2, lwd=4)

plotWAC(D$Time, D$GGVEWB, ylab="VEW")
lineWAC(D$Time, D$GGVNSB, col='darkgreen')
# define a chisquare error function:
wx <- 12.3
wy <- 13.2
V <- 154.
hdg <- D$THDG * Cradeg
csq <- function (x, D) {
  hdg <- (D$THDG+x[4]) * pi / 180.
  dvx <- D$GGVEWB - (x[2] + x[1]*sin(hdg))
  dvy <- D$GGVNSB - (x[3] + x[1]*cos(hdg))
  chisq <- sum (dvx**2 + dvy**2)
}
A <- nlm (csq, c(V, wx, wy, -0.1), D, hessian=TRUE)
print(A$estimate)
V <- A$estimate[1]
wx <- A$estimate[2]
wy <- A$estimate[3]
print(sprintf(" wind estimate: %.1f / %.1f", atan2(wx,wy)*180/pi+180, sqrt(wx**2+wy**2)))
vx <- V * cos(hdg) + wx
vy <- V * sin(hdg) + wy
lineWAC (D$Time, vx, col='red', lty=2, lwd=2)
lineWAC (D$Time, vy, col='orange', lty=2, lwd=2)
rms <- (A$minimum/length(D$THDG))**0.5

