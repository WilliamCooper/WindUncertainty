#<<initialization,echo=FALSE,include=FALSE>>=
require(Ranadu)
require(ggplot2)
require(grid)
library(knitr)
require(ggthemes)
require(vioplot)
require(plyr)
Directory <- "/home/Data/"      
Project = "DEEPWAVE"

# Flight <- "rf01"       
# fname = sprintf("%s%s/%s%s.nc", Directory,Project,Project,Flight)
# VarNames <- c("VYC","GGALT","LATC", "LONC", "PSXC", "QCXC",
#               "WDC", "WSC", "GGVEW", "GGVNS", "VEW", "VNS", "TASX",
#               "ADIFR", "AKRD", "SSLIP", "PITCH", "GGTRK", 
#               "ROLL", "THDG", "BDIFR", "EWX", "WIC",
#               "VSPD", "ATX", "GGVSPDB")
# Data <- getNetCDF (fname, VarNames, F=1) 
# for (i in c(3,4,6,7,8,9,10,11,12,13,15)) {
#   Flight <- sprintf("rf%02d",i)
#   fname = sprintf("/home/Data/%s/%s%s.nc", Project,Project,Flight)
#   Data <- merge(Data, getNetCDF (fname, VarNames, F=i), all=TRUE)
# }
# 
# Flight <- "rf01"       
# fname = sprintf("%s%s/%s%s.nc", Directory,Project,Project,Flight)
# VarNames <- c("VYC","GGALT", "PSXC", "QCXC", "WDC", "WSC", "GGVEW", 
#               "GGVNS", "VEW", "VNS", "TASX", "ADIFR", "AKRD", 
#               "SSLIP", "PITCH", "GGTRK", "ROLL", "THDG", "BDIFR", 
#               "EWX", "WIC", "VSPD", "ATX", "GGVSPDB", 
#               "VSPD_A", "CVSPD_GP", "CVSPD_GP")
# Data <- getNetCDF (fname, VarNames, F=1) 
# for (i in c(1:26)) {
#   Flight <- sprintf("rf%02d",i)
#   fname = sprintf("/home/Data/%s/%s%s.nc", Project,Project,Flight)
#   Data <- merge(Data, getNetCDF (fname, VarNames, F=i), all=TRUE)
# }
# save(Data,file=sprintf("%s%s/DWdata.Rdata.gz", Directory, Project), 
#      compress="gzip")
load(file=sprintf("%s%s/DWdata.Rdata.gz", Directory, Project))
  