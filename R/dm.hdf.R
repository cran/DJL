dm.hdf <-
function(xdata, ydata, rts = "crs", 
                   wd = NULL, se = FALSE, sg = "ssm", date = NULL, cv = "convex"){
  
  # Initial checks
  if(is.na(match(rts, c("crs", "vrs", "irs", "drs")))) stop('rts must be "crs", "vrs", "irs", or "drs".')
  if(is.na(match(se,  c(0, 1, FALSE, TRUE))))          stop('se must be either 0(FALSE) or 1(TRUE).')
  if(is.na(match(sg,  c("ssm", "max", "min"))))        stop('sg must be "ssm", "max", or "min".')
  if(is.na(match(cv,  c("convex", "fdh"))))            stop('cv must be "convex" or "fdh".')
  
  # Load library
  # library(lpSolveAPI)
  
  # Parameters
  xdata <- as.matrix(xdata)
  ydata <- as.matrix(ydata)
  g     <- cbind(xdata, ydata)
  date  <- if(!is.null(date)) as.matrix(date)
  n     <- nrow(xdata)
  m     <- ncol(xdata)
  s     <- ncol(ydata)
  wd    <- if(is.null(wd)) matrix(c(0), ncol = s) else as.matrix(wd)
  se    <- ifelse(is.logical(se), ifelse(isTRUE(se), 1, 0), se)
  rts   <- ifelse(cv == "fdh", "vrs", rts)
  
  # Data frames
  results.efficiency <- matrix(NA, nrow = n, ncol = 1)
  results.iteration  <- matrix(NA, nrow = n, ncol = 1) 
  
  # Hyperplane consistency checking function
  on_h <- function(g){
    r.h         <- dm.sf(xdata, ydata, rts, g, wd, se, sg, date, cv)
    w           <- r.h$w[k,]
    p           <- r.h$p[k,]
    u           <- r.h$u[k,]
    theta       <- ifelse(sum(w * xdata[k,]) == 0, -sum(p * ydata[k,]) / u, 
                          ifelse(sum(p * ydata[k,]) == 0, -u / sum(w * xdata[k,]),
                                 (-u + (u^2 - 4 * sum(w * xdata[k,]) * sum(p * ydata[k,]))^0.5) / (2 * sum(w * xdata[k,]))))
    z.xdata     <- xdata
    z.xdata[k,] <- xdata[k,] * theta
    z.ydata     <- ydata
    z.ydata[k,] <- ydata[k,] / theta
    z.g         <- cbind(z.xdata, z.ydata)
    on.h        <- isTRUE(round(abs(dm.sf(z.xdata, z.ydata, rts, z.g, wd, T, sg, date, cv)$eff[k]), 9) == 0)
    new.g       <- abs(c(xdata[k,] - theta * xdata[k,], ydata[k,]/theta - ydata[k,]))
    results     <- list(test = on.h, theta = theta, new.g = new.g)
    return(results)
  }
  
  # Run
  z.g <- g
  r   <- dm.sf(xdata, ydata, rts, z.g, wd, se, sg, date, cv)
  for (k in 1:n){
    if(r$eff[k] == 0){
      results.efficiency[k] <- 1
      results.iteration[k]  <- 1
    }else{
      h    <- 0
      flag <- T
      while(flag){
        z       <- on_h(z.g)
        flag    <- !z$test
        theta   <- z$theta
        z.g[k,] <- z$new.g
        h       <- h + 1
      }
      results.efficiency[k] <- theta
      results.iteration[k]  <- h
    }
  }
  temp.f <- dm.sf(xdata, ydata, rts, z.g, wd, se, sg, date, cv)
  
  # Get results
  results <- list(eff = results.efficiency, lambda = temp.f$lambda, mu = temp.f$mu,
                  xslack = temp.f$xslack, yslack = temp.f$yslack, iteration = results.iteration)
  return(results)
}
