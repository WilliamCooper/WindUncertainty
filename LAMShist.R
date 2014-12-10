require(Ranadu)
Flight <- "rf10"
SavePlotsToFiles <- FALSE
fname = sprintf("/home/Data/DEEPWAVE/DEEPWAVE%s.nc", Flight)
VarNames <- standardVariables()
VarNames <- c(VarNames, "ATHR1", "ATHR2", "ATRL", "AT_A", "AT_A2")
VarNames <- c(VarNames, "EW_DPL", "DP_DPL", "PALTF", "EW_DPR",
              "DP_DPR", "EW_VXL", "DP_VXL", "CAVP_CR2",
              "MIRRTMP_CR2")
VarNames <- c(VarNames, "PSFC", "PS_A", "PS_GP", "QCFC", 
              "QCRC", "QC_A", "QC_GP")
VarNames <- c(VarNames, "GGQUAL", "GGVEW", "GGVNS", "VEW", "VNS",
              "AKRD", "SSLIP", "THDG")
VarNames <- c(VarNames, "CONCU_RWO", "CONCU100_RWO", "CONCU500_RWO")
Data <- getNetCDF (fname, VarNames)
attach(Data)
netCDFfile = open.ncdf(fname)
V1 <- get.var.ncdf(netCDFfile, "BEAM1_LAMS")
V2 <- get.var.ncdf(netCDFfile, "BEAM2_LAMS")
V3 <- get.var.ncdf(netCDFfile, "BEAM3_LAMS")
close.ncdf (netCDFfile)

tck <-110000    #80904 looks reasonable; 80946 best rf03?
rg <- 1:512
for (k in seq(0,3000,100)) {
  t <- getIndex(Time, tck+k)
  plotWAC (rg, V1[rg,t], xlab="LAMS speeds", ylim=c(1.e7, 2.e9), log="y")
  points (rg, V2[rg,t], type='l', col="red")
  points (rg, V3[rg,t], type='l', col="darkgreen")
  legend ("topright", legend=c("Beam 1", "Beam 2", "Beam 3"), lty=1, 
        col=c("blue", "red", "darkgreen"), cex=0.8)
  title(sprintf("%.0f altitude %.0f ft TAS %.1f m/s", tck+k, GGALT[t]/0.3048, TASX[t]))
  print (GGALT[t]/0.3048)
}
if (SavePlotsToFiles) {
  dev.off()
}
detach (Data)
#rstudio::viewer("~/RStudio/DEEPWAVE/DWff03Plot03.png")




