Dc <- D[setRange (D$Time, startC, endC), ]
Dc <- Dc[abs(Dc$ROLL) > 26, ]
DR <- Dc[Dc$ROLL > 0, ]
DL <- Dc[Dc$ROLL < 0, ]
Aboth <- fitGS (csq, Dc)
rms <- (Aboth$minimum/nrow(Dc))**0.5
ALeft <- fitGS(csq, DL)
rmsL <- (ALeft$minimum/nrow(DL))**0.5
ARight <- fitGS (csq, DR)
rmsR <- (ARight$minimum/nrow(DR))**0.5
bestFit <- Aboth$estimate
bestFitL <- ALeft$estimate
bestFitR <- ARight$estimate
bestWD <- (atan2(bestFit[2], bestFit[3]) / Cradeg) %% 360
bestWS <- sqrt(bestFit[2]**2 + bestFit[3]**2)
mwd <- mean(Dc$WDC, na.rm=TRUE)
mws <- mean(Dc$WSC, na.rm=TRUE)
mtas <- mean(Dc$TASX, na.rm=TRUE)
bestWDL <- (atan2(bestFitL[2], bestFitL[3]) / Cradeg) %% 360
bestWSL <- sqrt(bestFitL[2]**2 + bestFitL[3]**2)
mwdL <- mean(DL$WDC, na.rm=TRUE)
mwsL <- mean(DL$WSC, na.rm=TRUE)
mtasL <- mean(DL$TASX, na.rm=TRUE)
bestWDR <- (atan2(bestFitR[2], bestFitR[3]) / Cradeg) %% 360
bestWSR <- sqrt(bestFitR[2]**2 + bestFitR[3]**2)
mwdR <- mean(DR$WDC, na.rm=TRUE)
mwsR <- mean(DR$WSC, na.rm=TRUE)
mtasR <- mean(DR$TASX, na.rm=TRUE)
