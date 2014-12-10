require(Ranadu)
require (ggplot2)
require(grid)
Flight <- "rf15"
Directory <- "/home/Data/"
Project <- "DEEPWAVE"
fname = sprintf("%s%s/DEEPWAVE%s.nc", Directory, Project, Flight)
VarNames <- c("VYC","GGALT","LATC", "LONC", "PSXC", "QCXC",
              "WDC", "WSC", "GGVEW", "GGVNS", "VEW", "VNS",
              "ADIFR", "AKRD", "SSLIP", "PITCH", "TASX",
              "ROLL", "THDG", "BDIFR", "EWX",
              "ADIFR_GP", "BDIFR_GP", "PS_GP", "QC_GP",
              "CROLL_GP", "CPITCH_GP", "CTHDG_GP", "WIC",
              "CVNS_GP", "CVEW_GP", "VSPD", "CVSPD_GP",
              "ATX")
Cradeg <- pi / 180.
D <- getNetCDF (fname, VarNames)      # this handles high-rate data also
Flight <- "rf11"        # this had  cal maneuvers at 40K ft
fname = sprintf("%s%s/DEEPWAVE%s.nc", Directory, Project, Flight)
D2 <- getNetCDF (fname, VarNames)
r4 <- setRange (D2$Time, 103000,104000)
D2 <- D2[r4,]
Flight <- "rf16"        # 43000-ft leg
fname = sprintf("%s%s/DEEPWAVE%s.nc", Directory, Project, Flight)
D3 <- getNetCDF (fname, VarNames)
r5 <- setRange (D3$Time, 94000,110000)
D3 <- D3[r5,]
Flight <- "rf14"        # 45000-ft leg
fname = sprintf("%s%s/DEEPWAVE%s.nc", Directory, Project, Flight)
D4 <- getNetCDF (fname, VarNames)
r6 <- setRange (D4$Time, 113000,122000)
D4 <- D4[r6,]
flight = "rf11"
fname = sprintf("%s%s/DEEPWAVE%s.nc", Directory, Project, Flight)
D5 <- getNetCDF (fname, VarNames)
r7 <- setRange (D5$Time, 70000, 100000)
D5 <- D5[r7,]
r1 <- setRange (D$Time, 32100,32900)
r2 <- setRange (D$Time, 41500,42300)
r3 <- setRange (D$Time, 50100,51100)
# need to add high-altitude segment also to constrain fit:
# rf11 103000--104000 is a speed run at 40,000 ft.
# useful also to include some data from max altitude:
r <- c(r1,r2,r3)
DD <- merge(D2, D4, all=T)
DD <- merge(DD, D[r,], all=TRUE)
DD <- merge (DD, D3, all=TRUE)
Data <- merge (DD, D5, all=TRUE)
Data2 <- merge (D2, D[r,], all=T)
Valid <- (Data2$TASX > 130.)
Data2 <- Data2[Valid,]
attach(Data2)
Cradeg <- pi / 180.
Mach <- MachNumber (PSXC, QCXC)  # uses conventional q, p
AOAREF <- PITCH - asin(VSPD/TASX) / Cradeg
AQR <- ADIFR/QCXC # basic pressure ratio for AOA
fmy <- lm(AOAREF~AQR)
summary(fmy)
cfr <- coefficients (fmy)
A1 <- cfr[1]+cfr[2]*AQR
plot (AOAREF, A1, pch=16, cex=0.8, col='blue', xlab="Ref. AOA", ylab="fit AOA")
lines (c(-3.,6.), c(-3.,6.), lty=2, lwd=3, col='darkorange')
AOAREF_GP <- CPITCH_GP - asin(CVSPD_GP/TASX) / Cradeg
AQR_GP <- ADIFR_GP/QC_GP # basic pressure ratio for AOA
RR2_GP <- QC_GP/PS_GP    # q/p ratio from the gust pod
#fm <- lm (AOAREF~AQR+RR2+PS_GP+QC_GP+Mach)
#fm_GP <- lm (AOAREF_GP~AQR_GP+RR2_GP+Mach)
AQRM_GP <- AQR_GP * Mach
fm_GP <- lm (AOAREF_GP~AQR_GP+AQRM_GP +RR2_GP)
cf3_GP <- coefficients (fm_GP)
#print(c("Fit coefficients are: ",cf3_GP))
print(summary(fm_GP))
A2_GP <- cf3_GP[1]+cf3_GP[2]*AQR_GP+cf3_GP[3]*AQRM_GP+cf3_GP[4]*RR2_GP
plot (AOAREF_GP, A2_GP, xlab="Predicted AOA [deg.]", 
      ylab="Fit AOA [deg.]", type='p', pch=16, cex=0.7, col='blue')
lines (c(-3.,6.), c(-3.,6.), lty=2, lwd=3, col='orange')
title ("Gust Pod")
theme_WAC <- theme_gdocs() + theme(
  plot.title = element_text(hjust=0.5),
  axis.text = element_text(size = 16),
  panel.grid.major = element_line(color = "lightblue", linetype=5),
  panel.background = element_rect(fill = "gray95"),
  axis.title=element_text(face="plain",size=18,colour="blue"),
  line=element_line(size=1),
  axis.ticks = element_line(size=1),
  axis.ticks.length = unit(-0.35,"cm"),
  axis.ticks.margin = unit (0.6,"cm"),
  legend.position=c(0.5,0.96),
  plot.margin=unit(c(1.5,1,0.5,0.5),"lines"),
  plot.title=element_text(vjust=1.3),
  legend.background=element_rect(fill="white"),
  legend.direction="horizontal",
  panel.border=element_rect(colour="black",size=0.7),
  axis.title.y=element_text(angle=90))
xp <- AOAREF_GP
yp <- A2_GP
g <- ggplot(data=Data2,aes(x=xp))
g <- g + geom_point(aes(y=yp), colour="blue", pch=19)
g <- g + xlab ("Predicted AOA [deg.]")
g <- g + ylab ("AOA from fit [deg.]") 
#g <- g + geom_line(aes(y=xp), colour="darkorange", lwd=1, linetype=2)
g <- g + geom_smooth(aes(y=yp),method="lm", colour="darkorange", lwd=2, linetype=2)
g <- g + theme_WAC
# g <- g + annotate("rect", xmin = -2.5, xmax = 0, ymin = 1, ymax = 2,
#              alpha = 1, fill="lightyellow")
# lbl <- expression(paste(Delta,p[alpha]))
# lbl <- expression(Delta)
# lbl <- c("Test\nTest2")
# lbl <- expression(alpha)
# g <- g + annotate("text", x = -2, y = 1.8, label = lbl)
g

# sideslip
detach(Data2)
# first get filtered wind-component values:
# u <- -1. * D$WSC * sin (D$WDC*Cradeg)
# v <- -1. * D$WSC * cos (D$WDC*Cradeg)
# u <- ButterworthFilter (u, tau=30)
# v <- ButterworthFilter (v, tau=30)

r1 <- setRange (D$Time, 33230,33500)
r2 <- setRange (D$Time, 43100,43300)
r3 <- setRange (D$Time, 52640,52930)
r <- c(r1,r2,r3)
# u <- u[r]
# v <- v[r]
Data3 <- D[r,]
Valid <- (Data3$TASX > 130.)
Data3 <- Data3[Valid,]
attach(Data3)
u <- -1. * WSC * sin (WDC*Cradeg)
v <- -1. * WSC * cos (WDC*Cradeg)
SSREF_GP <- -CTHDG_GP + atan2((CVEW_GP-u), (CVNS_GP-v))/ Cradeg
SSREF_GP[SSREF_GP < -180.] <- SSREF_GP[SSREF_GP < -180.] + 360.
BQR_GP <- BDIFR_GP / QC_GP
Mach <- MachNumber (PSXC, QCXC)
RR2_GP <- QC_GP / PS_GP
sfm_GP <- lm(SSREF_GP~BQR_GP)
summary(sfm_GP)
cfs_GP <- coefficients(sfm_GP)
S2_GP <- cfs_GP[1] + cfs_GP[2] * BQR_GP 
plot (SSREF_GP, S2_GP, pch=16, cex=0.6, col='blue',
      xlab="SS Reference", ylab="SSLIP from fit")
lines(c(-2.,2.), c(-2.,2.), col='orange', lty=2, lwd=3)
title ("Gust Pod")
detach(Data3)
xp <- SSREF_GP
yp <- S2_GP
g <- ggplot(data=Data3,aes(x=xp))
g <- g + geom_point(aes(y=yp), colour="blue", pch=19)
g <- g + xlab ("Predicted sideslip [deg.]")
g <- g + ylab ("Sideslip from fit [deg.]") 
#g <- g + geom_line(aes(y=xp), colour="darkorange", lwd=1, linetype=2)
g <- g + geom_smooth(aes(y=yp),method="lm", colour="darkorange", lwd=1.0, linetype=2)
g <- g + theme_WAC
g

#TAS
# for this purpose, use the dataset including speed runs and other-altitude segments
#Data <- D
Valid <- (Data$TASX > 130.)
Data <- Data[Valid,]
attach(Data)
B1 <- QCXC/PSXC
B2 <- QC_GP/PS_GP
AQR_GP <- ADIFR_GP/QC_GP # basic pressure ratio for AOA
RR2_GP <- QC_GP/PS_GP    # q/p ratio from the gust pod
Mach <- MachNumber(PS_GP,QC_GP)
# B2 and AQR are based only on the gust-pod measurements
#fm4 <- lm(B1~B2+I(B2^2)+AQR_GP)
BxM <- B2*Mach
fm4 <- lm(B1~B2+BxM+I(B2^2)+AQR_GP+Mach)
cf4 <- coefficients (fm4)
#print (c("fit coefficients for TAS:", cf4))
print(summary(fm4))
#B3 <- cf4[1]+cf4[2]*B2+cf4[3]*B2**2+cf4[4]*AQR_GP
B3 <- cf4[1]+cf4[2]*B2+cf4[3]*BxM+cf4[4]*B2**2+cf4[5]*AQR_GP+cf4[6]*Mach

# B3 is the deduced pressure ratio to use when calculating TAS
plot (B1, B3, xlab="Standard q/p", ylab="gust-probe q/p", xlim=c(0.1,0.6), ylim=c(0.1,0.6), pch=16, cex=0.7, col='blue')
lines (c(0., 1.), c(0.,1.), lwd=3, lty=2, col='orange')
xp <- B1
yp <- B3
g <- ggplot(data=Data,aes(x=xp))
g <- g + geom_point(aes(y=yp), colour="blue", pch=19)
g <- g + xlab ("Standard q/p")
g <- g + ylab ("Fit q/p from gust-pod variables") 
#g <- g + geom_line(aes(y=xp), colour="darkorange", lwd=1, linetype=2)
g <- g + geom_smooth(aes(y=yp),method="lm", colour="darkorange", lwd=1.0, linetype=2)
g <- g + theme_WAC
g

Vg <- TrueAirspeed (MachNumber (PS_GP, PS_GP*B3, EWX), ATX)
plot(TASX, Vg, xlab="TASX [m/s]", ylab="TAS, gust pod [m/s]",
     xlim=c(200.,240.), ylim=c(200.,240.),pch=16, cex=0.7, 
     col='blue')
lines (c(100.,300.), c(100.,300.), lty=2, col='orange')
xp <- TASX
yp <- Vg
g <- ggplot(data=Data,aes(x=xp))
g <- g + geom_point(aes(y=yp), colour="blue", pch=19)
g <- g + xlab ("TASX [m/s]")
g <- g + ylab ("TASG [m/s]") 
#g <- g + geom_line(aes(y=xp), colour="darkorange", lwd=1, linetype=2)
g <- g + geom_smooth(aes(y=yp),method="lm", colour="darkorange", lwd=2.0, linetype=2)
g <- g + theme_WAC
g
detach(Data)
# processor:
# process a specific flight:
Flight <- "rf16"
Directory <- "/home/Data/"
Project <- "DEEPWAVE"
fname = sprintf("%s%s/DEEPWAVE%s.nc", Directory, Project,Flight)
# copy to a new file before adding variables to it:
fnew = sprintf("%s%s/DW%s_GP.nc", Directory, Project,Flight)
#system(paste ("cp", fname, fnew, sep=' '), wait=TRUE)
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
Valid <- (Data$TASX > 130.)
Valid <- (Valid & ((abs(Data$ROLL) < 2.) & (abs(Data$VSPD) < 4.)))
DataV <- Data[Valid,]  # get the valid-points data file
attach(Data)           # attach the full data file
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
# instead, use TASX:
TAS_GP <- TASX
# define a dataframe for the relative wind from the gust pod:
# d <- data.frame("U_RW"=TAS_GP)
# d["V_RW"] <- TAS_GP * tan (SS_GP * Cradeg)
# d["W_RW"] <- TAS_GP * tan (AOA_GP * Cradeg)
TAS <- TASX / sqrt(1.+tan(AOA_GP*Cradeg)**2+tan(SS_GP*Cradeg))
d <- data.frame("U_RW"=TAS)
d["V_RW"] <- TAS * tan (SS_GP * Cradeg)
d["W_RW"] <- TAS * tan (AOA_GP * Cradeg)
TAS <- TASX / sqrt(1.+tan(AOA*Cradeg)**2+tan(SS*Cradeg))
d2 <- data.frame ("U2"=TAS)
d2["V2"] <- TAS * tan (SS * Cradeg)
d2["W2"] <- TAS * tan (AOA * Cradeg)
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
T1 <- aperm(array (c(One,Z,Z,Z,cosphi,-sinphi,Z,sinphi,cosphi), 
                   dim=c(DL,3,3)))
T2 <- aperm(array (c (costheta,Z,sintheta,Z,One,Z,-sintheta,Z,costheta), 
                   dim=c(DL,3,3)))
T3 <- aperm(array (c (cospsi,-sinpsi,Z,sinpsi,cospsi,Z,Z,Z,One), 
                   dim=c(DL,3,3)))
TT1 <- aperm(array (c(One,Z,Z,Z,cosphi2,-sinphi2,Z,sinphi2,cosphi2), 
                    dim=c(DL,3,3)))
TT2 <- aperm(array (c (costheta2,Z,sintheta2,Z,One,Z,-sintheta2,Z,costheta2), 
                    dim=c(DL,3,3)))
TT3 <- aperm(array (c (cospsi2,-sinpsi2,Z,sinpsi2,cospsi2,Z,Z,Z,One), 
                    dim=c(DL,3,3)))
WDG <- vector ("numeric", DL)
WSG <- vector ("numeric", DL)
WIG <- vector ("numeric", DL)
WDX <- vector ("numeric", DL)
WSX <- vector ("numeric", DL)
WIX <- vector ("numeric", DL)
TASG <- vector ("numeric", DL)
CVEWC_GP <- ComplementaryFilter (CVEW_GP, GGVEW, 150.)
CVNSC_GP <- ComplementaryFilter (CVNS_GP, GGVNS, 150.)
Hlast <- 0.

# here is the loop equivalent:
for (i in 1:DL) {
  Y1 <- T1[,,i] %*% matrix (rw[i,], 3, 1)
  Y2 <- T2[,,i] %*% Y1
  Y3 <- T3[,,i] %*% Y2
  WG <- matrix (c(-CVNSC_GP[i], -CVEWC_GP[i], CVSPD_GP[i]), 
                3, 1)
  Y4 <- Y3 + WG
  WDG[i] <- atan2 (Y4[2], Y4[1]) / Cradeg
  if ((!is.na(WDG[i])) & (WDG[i] < 0.)) {
    WDG[i] <- WDG[i] + 360.
  }
  WSG[i] <- sqrt (Y4[1]**2 + Y4[2]**2)
  WIG[i] <- Y4[3]
  # attempt to flag bad CTHDG_GP points as missing
  if (!is.na(CTHDG_GP[i])) {
    if ((abs(CTHDG_GP[i]-Hlast) > 3.) 
        & (abs(Hlast-180.) < 5.)) {
      WIG[i] <- NA
      WDG[i] <- NA
      WSG[i] <- NA
    }
  } else {
    WIG[i] <- NA
    WDG[i] <- NA
    WSG[i] <- NA    
  }
  if (!is.na(CTHDG_GP[i])) {
    Hlast <- CTHDG_GP[i]
  }
  Y1 <- TT1[,,i] %*% matrix (rw2[i,], 3, 1)
  Y2 <- TT2[,,i] %*% Y1
  Y3 <- TT3[,,i] %*% Y2
  WG <- matrix (c(-VNS[i], -VEW[i], VSPD[i]), 
                3, 1)
  Y4 <- Y3 + WG
  WDX[i] <- atan2 (Y4[2], Y4[1]) / Cradeg
  if ((!is.na(WDX[i])) & (WDX[i] < 0.)) {
    WDX[i] <- WDX[i] + 360.
  }
  WSX[i] <- sqrt (Y4[1]**2 + Y4[2]**2)
  WIX[i] <- Y4[3]
}

detach(Data)
Data["WIX"] <- WIX
Data["WIG"] <- WIG
Data["WDX"] <- WDX
Data["WSX"] <- WSX
Data["WDG"] <- WDG
Data["WSG"] <- WSG
Data["AOA_GP"] <- AOA_GP
Data["SS_GP"] <- SS_GP
rm(WIX,WIG,WDX,WDG,WSX,WSG,AOA_GP,SS_GP)
DataV <- Data[Valid,]

attach(DataV)
xp <- WIX
yp <- WIG
g <- ggplot(data=DataV,aes(x=WIX))
g <- g + geom_point(aes(y=WIG), colour="blue", pch=19)
g <- g + xlab ("WIX [m/s]")
g <- g + ylab ("WIG [m/s]") 
#g <- g + geom_line(aes(y=xp), colour="darkorange", lwd=1, linetype=2)
g <- g + geom_smooth(aes(y=yp),method="lm", colour="darkorange", lwd=2.0, linetype=2)
g <- g + theme_WAC
g
print(summary(lm(WIG~WIX)))
xp <- AKRD
yp <- AOA_GP
g <- ggplot(data=DataV,aes(x=AKRD))
g <- g + geom_point(aes(y=AOA_GP), colour="blue", pch=19)
g <- g + xlab ("AKRD [deg.]")
g <- g + ylab ("AOA from Gust Pod [deg.]") 
#g <- g + geom_line(aes(y=xp), colour="darkorange", lwd=1, linetype=2)
#g <- g + geom_smooth(aes(y=yp),method="lm", colour="darkorange", lwd=2.0, linetype=2)
g <- g + theme_WAC
g
print(summary(lm(AOA_GP~AKRD)))
print(summary(lm(AKRD~AOA_GP)))
xpr <- AKRD-mean(AKRD, na.rm=TRUE)
ypr <- AOA_GP-mean(AOA_GP, na.rm=TRUE)
xpr2 <- mean(xpr**2,na.rm=TRUE)
ypr2 <- mean(ypr**2,na.rm=TRUE)
xprypr <- mean(xpr*ypr, na.rm=TRUE)
tan2theta <- 2.*xprypr / (xpr2-ypr2)
twotheta <- atan(tan2theta)
b <- tan(twotheta/2.)
if (b*xprypr < 0.) {
  b <- -1./tan(twotheta/2.)
}
yzero <- mean(AOA_GP, na.rm=TRUE)-b*mean(AKRD, na.rm=TRUE)
zp <- yzero + b*xp
print (c(yzero,b))
g + geom_line(aes(y=zp),colour="darkorange", lwd=2)
detach(DataV)
attach(Data)

clr <- c("WDX","WDG")
col <- c("red","blue")
op <- par (mfrow=c(2,1), mar=c(4,4,0,2)+0.1)
xp <- Time
r <- setRange(Time, 60000,125000)
#r <- setRange(Time,0,0)
DataP <- Data[r,]
g <- ggplot(data=DataP,aes(x=Time))
g <- g + geom_line(aes(y=WDX, colour="WDX"), lwd=1)
g <- g + geom_line (aes(y=WDG, colour="WDG"), lwd=1)
g <- g + xlab ("Time [UTC]")
g <- g + ylab ("Wind Direction [deg.]") 
#g <- g + geom_line(aes(y=xp), colour="darkorange", lwd=1, linetype=2)
#g <- g + geom_smooth(aes(y=yp),method="lm", colour="darkorange", lwd=2.0, linetype=2)
g <- g + scale_colour_manual("", 
                              breaks = clr,
                              values = col)
clr <- c("WSX","WSG")
g1 <- g + ylim(200.,300.) + theme_WAC
g <- ggplot(data=DataP,aes(x=Time))
g <- g + geom_line(aes(y=WSX, colour="WSX"), lwd=1)
g <- g + geom_line (aes(y=WSG, colour="WSG"), lwd=1)
g <- g + xlab ("Time [UTC]")
g <- g + ylab ("Wind Speed [m/s]") 
#g <- g + geom_line(aes(y=xp), colour="darkorange", lwd=1, linetype=2)
#g <- g + geom_smooth(aes(y=yp),method="lm", colour="darkorange", lwd=2.0, linetype=2)
g <- g + scale_colour_manual("", 
                             breaks = clr,
                             values = col)
g2 <- g + theme_WAC
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(g1,g2,cols=1)

clr <- c("WIX","WIG")
r <- setRange(Time, 83000,85000)
#r <- setRange(Time,0,0)
DataP <- Data[r,]
g <- ggplot(data=DataP,aes(x=Time))
g <- g + geom_line(aes(y=WIX, colour="WIX"), lwd=1.2)
g <- g + geom_line (aes(y=WIG, colour="WIG"), lwd=0.7)
g <- g + xlab ("Time [UTC]")
g <- g + ylab ("Vertical Wind [m/s]") 
#g <- g + geom_line(aes(y=xp), colour="darkorange", lwd=1, linetype=2)
#g <- g + geom_smooth(aes(y=yp),method="lm", colour="darkorange", lwd=2.0, linetype=2)
g <- g + scale_colour_manual("", 
                             breaks = clr,
                             values = col)
g <- g + theme_WAC
g

clr <- c("ROLL","CROLL_GP")
col <- c("red","blue")
xp <- Time
r <- setRange(Time, 84900,85820)
#r <- setRange(Time,0,0)
DataP <- Data[r,]
g <- ggplot(data=DataP,aes(x=Time))
g <- g + geom_line(aes(y=ROLL, colour="ROLL"), lwd=1)
g <- g + geom_line (aes(y=CROLL_GP, colour="CROLL_GP"), lwd=1)
g <- g + xlab ("Time [UTC]")
g <- g + ylab ("Roll [deg.]") 
#g <- g + geom_line(aes(y=xp), colour="darkorange", lwd=1, linetype=2)
#g <- g + geom_smooth(aes(y=yp),method="lm", colour="darkorange", lwd=2.0, linetype=2)
g1 <- g + scale_colour_manual("", 
                             breaks = clr,
                             values = col) + theme_WAC + theme(legend.position=c(0.5,1.15))
clr <- c("THDG","CTHDG_GP")
col <- c("red","blue")
xp <- Time
g <- ggplot(data=DataP,aes(x=Time))
g <- g + geom_line(aes(y=THDG, colour="THDG"), lwd=1)
g <- g + geom_line (aes(y=CTHDG_GP, colour="CTHDG_GP"), lwd=1)
g <- g + xlab ("Time [UTC]")
g <- g + ylab ("Heading [deg.]") 
#g <- g + geom_line(aes(y=xp), colour="darkorange", lwd=1, linetype=2)
#g <- g + geom_smooth(aes(y=yp),method="lm", colour="darkorange", lwd=2.0, linetype=2)
g2 <- g + scale_colour_manual("", 
                              breaks = clr,
                              values = col) + theme_WAC
clr <- c("PITCH","CPITCH_GP")
col <- c("red","blue")
xp <- Time
g <- ggplot(data=DataP,aes(x=Time))
g <- g + geom_line(aes(y=PITCH, colour="PITCH"), lwd=1)
g <- g + geom_line (aes(y=CPITCH_GP, colour="CPITCH_GP"), lwd=1)
g <- g + xlab ("Time [UTC]")
g <- g + ylab ("Pitch [deg.]") 
#g <- g + geom_line(aes(y=xp), colour="darkorange", lwd=1, linetype=2)
#g <- g + geom_smooth(aes(y=yp),method="lm", colour="darkorange", lwd=2.0, linetype=2)
g3 <- g + scale_colour_manual("", 
                              breaks = clr,
                              values = col) + theme_WAC + theme(legend.position=c(0.5,1.15))
multiplot(g1,g2,g3,cols=1)

detach(Data)


