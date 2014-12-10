require(Ranadu)
require(ggplot2)
require(grid)
library(knitr)
require(ggthemes)
require(vioplot)
require(plyr)
Directory <- DataDirectory()
Project = "DEEPWAVE"
VarNames <- c("TASX", "PITCH", "THDG", "GGALT", "GGVSPDB", "ADIFR", "QCF", "CPITCH_GP", "CPITCH_LAMS",
              "PSF", "QCF", "PITCH_IRS2", "PITCH_IRS3")
# assemble a data file for DEEPWAVE to use for the calibration
Flight <- "rf01"
NumberOfFlights <- 26
fname = sprintf("%s%s/%s%s.nc", Directory,Project,Project,Flight)
# Data <- getNetCDF (fname, VarNames, F=1) 
# for (i in c(2:NumberOfFlights)) {
#   Flight <- sprintf("rf%02d",i)
#   fname = sprintf("%s%s/%s%s.nc", Directory,Project,Project,Flight)
#   Data <- merge(Data, getNetCDF (fname, VarNames, F=i), all=TRUE)
# }
# save(Data,file="~/RStudio/DEEPWAVE/DWpitch.Rdata.gz", compress="gzip")
load(file="~/RStudio/DEEPWAVE/DWpitch.Rdata.gz")
# exclude obviously bad periods:
Cradeg = pi / 180
Valid <- Data$TASX > 130.
Valid[(Data$RF == 5)][setRange(Data$Time[Data$RF == 5], 0, 80000)] <- FALSE
Valid[Data$RF == 6] <- FALSE
Valid[Data$RF == 7] <- FALSE
Valid[(Data$RF == 8)][setRange(Data$Time[Data$RF == 8], 130000, 0)] <- FALSE
Valid[(Data$RF == 10)][setRange(Data$Time[Data$RF == 10], 0, 61500)] <- FALSE
Valid[(Data$RF == 15)][setRange(Data$Time[Data$RF == 15], 0, 30000)] <- FALSE
Valid[(Data$RF == 17)][setRange(Data$Time[Data$RF == 17], 0, 81000)] <- FALSE
Valid[(Data$RF == 20)][setRange(Data$Time[Data$RF == 20], 41500, 52000)] <- FALSE
Valid[(Data$RF == 23)][setRange(Data$Time[Data$RF == 23], 100000,103000)] <- FALSE
Valid[(Data$RF == 25)][setRange(Data$Time[Data$RF == 25], 0,63000)] <- FALSE
Valid[abs(Data$ROLL) > 5] <- FALSE
DV <- Data[Valid,]
DV["Mach"] <- MachNumber (DV$PSF, DV$QCF)  # uses uncorrected q, p
DV["AOAREF"] <- DV$PITCH - asin(DV$GGVSPDB/DV$TASX) / Cradeg
DV["AQR"] <- DV$ADIFR/DV$QCF # basic pressure ratio for AOA
DV["AQRM"] <- DV$AQR * DV$Mach
DV["RR2"] <- DV$QCF/DV$PSF    # q/p ratio
# for (i in 1:NumberOfFlights) {
#   if (i == 6) {next}
#   if (i == 7) {next}
#   t <- DV$RF == i
#   plot (DV$PITCH[t], DV$CPITCH_GP[t], pch=20, col='blue')
#   points (DV$PITCH[t], DV$CPITCH_LAMS[t], pch=20, col='red')
#   title (sprintf("Flight %d", i))
# }
t <- DV$RF == 12
D <- DV[t, ]
D <- na.omit(D)
c0 <- 4.45543
c1 <- 21.40075
D["AKRD"] <- c0 + c1 * (D$ADIFR/D$QCF)
D["WIX"] <- D$TASX * (D$AKRD-D$PITCH) * Cradeg + D$GGVSPDB
ww <- filter(butter(3,0.5/500), D$WIX)
plotWAC(D$Time, D$WIX, ylab='WIX')
lineWAC(D$Time, ww, col='red')
