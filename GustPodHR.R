require(Ranadu)

Flight <- "rf05"
SavePlotsToFiles <- TRUE
fname = sprintf("/Data/DEEPWAVE/DW%shr_GP.nc", Flight)

VarNames <- c("VYC","GGALT","LATC", "LONC", "PSXC", "QCXC",
              "WDC", "WSC", "GGVEW", "GGVNS", "VEW", "VNS", "TASX",
              "AKRD", "SSLIP", "PITCH", "ROLL", "THDG",
              "ADIFR_GP", "BDIFR_GP", "PS_GP", "QC_GP",
              "CROLL_GP", "CPITCH_GP", "CTHDG_GP", "WIC",
              "CVNS_GP", "CVEW_GP", "VSPD", "CVSPD_GP",
              "ATX", "WDG", "WSG", "WIG")
#VarNames <- standardVariables(vn)

D <- getNetCDF (fname, VarNames)
r <- setRange (D$Time, 91500,94000)
Data <- D[r,]
attach(Data)
plotWAC (Time, WIG)
# Best match is with a WIG shift earlier by 1/25-s:
N <-2
WIG2 <- c(WIG[N:length(WIG)], rep(0.,N-1))
plot (WIC, WIG2, pch=16, cex=0.6, col='blue')
r <- 1:(length(WIG)-N+1)
fm <- lm(WIG2[r]~WIC[r])
summary(fm)
WIG <- WIG2
plotWAC (Time, WIG)
points (Time, WIC, type='l', col='red', lwd=2)
detach(Data)

