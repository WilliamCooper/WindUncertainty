## @knitr getGGshift
getGGshift <- function (D, csqf) {  
  ## D=data.frame with GGVNS, GGVEW and THDG
  ## csq=minimization function for fitting
  NL <- dim(D)[1]
  rmsmin <- 1e6
  steps <- 30
  nbest <- 0
  CRadeg <- pi/180
  DSGGVNS <- D$GGVNS
  DSGGVEW <- D$GGVEW
  for (n in 0:steps) {
    D$GGVNS <- c(D$GGVNS[(1+n):NL],rep(D$GGVNS[NL],n))
    D$GGVEW <- c(D$GGVEW[(1+n):NL],rep(D$GGVEW[NL],n))
    wx <- 12.3
    wy <- 13.2
    V <- 154.
    hdg <- D$THDG * Cradeg
    A <- nlm (csqf, c(V, wx, wy, -0.1), D, hessian=TRUE)
    D$GGVEW <- DSGGVEW
    D$GGVNS <- DSGGVNS
    #print(A$estimate)
    V <- A$estimate[1]
    wx <- A$estimate[2]
    wy <- A$estimate[3]
    dh <- A$estimate[4]
    #print(sprintf(" wind estimate: %.1f / %.1f", atan2(wy,wx)*180/pi+180, sqrt(wx**2+wy**2)))
    vx <- V * cos(hdg+dh*Cradeg) + wx
    vy <- V * sin(hdg+dh*Cradeg) + wy
    rms <- (A$minimum/length(D$THDG))**0.5
    ##print(sprintf("RMS error for shift of %d samples: %5.2f m/s", n, rms))
    if (rms < rmsmin) {
      rmsmin <- rms
      bestFit <- A$estimate
      nbest <- n
    }
  }
  for (n in 1:steps) {    # then shift backward
    D$GGVNS <- c(rep(D$GGVNS[1],n), D$GGVNS[1:(NL-n)])
    D$GGVEW <- c(rep(D$GGVEW[1],n), D$GGVEW[1:(NL-n)])
    wx <- 12.3
    wy <- 13.2
    V <- 154.
    hdg <- D$THDG * Cradeg
    A <- nlm (csqf, c(V, wx, wy, -0.1), D, hessian=TRUE)
    D$GGVEW <- DSGGVEW
    D$GGVNS <- DSGGVNS
    #print(A$estimate)
    V <- A$estimate[1]
    wx <- A$estimate[2]
    wy <- A$estimate[3]
    dh <- A$estimate[4]
    #print(sprintf(" wind estimate: %.1f / %.1f", atan2(wy,wx)*180/pi+180, sqrt(wx**2+wy**2)))
    vx <- V * cos(hdg+dh*Cradeg) + wx
    vy <- V * sin(hdg+dh*Cradeg) + wy
    rms <- (A$minimum/length(D$THDG))**0.5
    if (rms < rmsmin) {
      rmsmin <- rms
      bestFit <- A$estimate
      nbest <- -n
    }
    ##print(sprintf("RMS error for shift of %d samples: %5.2f m/s", -n, rms))
  }
  ##print (sprintf ("best shift is %d samples", nbest))
  return(c(nbest, rmsmin, bestFit))
}

