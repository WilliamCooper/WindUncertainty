Dc$xi <- Cradeg * (((Dc$THDG + Dc$SSLIP * cos (Dc$ROLL*Cradeg) - Dc$AKRD * sin (Dc$ROLL*Cradeg)-mwd)+360) %% 360)
fmcircleC <- lm (WSC ~ I(cos(xi))+I(TASX*sin(xi)), data=Dc)
# print (summary (fmcircle))
cfcircleC <- coefficients (fmcircleC)
VerrorC <- -cfcircleC[2]
HerrorC <- -cfcircleC[3] / Cradeg
DL$xi <- Cradeg * (((DL$THDG + DL$SSLIP * cos (DL$ROLL*Cradeg) - DL$AKRD * sin (DL$ROLL*Cradeg)-mwdL)+360) %% 360)
fmcircleL <- lm (WSC ~ I(cos(xi))+I(TASX*sin(xi)), data=DL)
# print (summary (fmcircle))
cfcircleL <- coefficients (fmcircleL)
VerrorL <- -cfcircleL[2]
HerrorL <- -cfcircleL[3] / Cradeg
DR$xi <- Cradeg * (((DR$THDG + DR$SSLIP * cos (DR$ROLL*Cradeg) - DR$AKRD * sin (DR$ROLL*Cradeg)-mwdR)+360) %% 360)
fmcircleR <- lm (WSC ~ I(cos(xi))+I(TASX*sin(xi)), data=DR)
# print (summary (fmcircle))
cfcircleR <- coefficients (fmcircleR)
VerrorR <- -cfcircleR[2]
HerrorR <- -cfcircleR[3] / Cradeg
