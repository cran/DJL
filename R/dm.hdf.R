dm.hdf <-
function(xdata, ydata, rts = "crs", 
                   wd = NULL, se = FALSE, sg = "ssm", date = NULL, cv = "convex", o = NULL){
  
  # Initial checks
  if(is.na(match(rts, c("crs", "vrs", "irs", "drs")))) stop('rts must be "crs", "vrs", "irs", or "drs".')
  if(is.na(match(se,  c(0, 1, FALSE, TRUE))))          stop('se must be either 0(FALSE) or 1(TRUE).')
  if(is.na(match(sg,  c("ssm", "max", "min"))))        stop('sg must be "ssm", "max", or "min".')
  if(is.na(match(cv,  c("convex", "fdh"))))            stop('cv must be "convex" or "fdh".')
  if(!is.null(o) && !all(o <= nrow(xdata)))            stop('o must be element(s) of n.')
  
  # Load library
  # library(lpSolveAPI)
  
  # Parameters
  xdata <- as.matrix(xdata)
  ydata <- as.matrix(ydata)
  date  <- if(!is.null(date)) as.matrix(date)
  n     <- nrow(xdata)
  m     <- ncol(xdata)
  s     <- ncol(ydata)
  wd    <- if(is.null(wd)) matrix(c(0), ncol = s) else as.matrix(wd)
  se    <- ifelse(is.logical(se), ifelse(isTRUE(se), 1, 0), se)
  rts   <- ifelse(cv == "fdh", "vrs", rts)
  o     <- if(is.null(o)) c(1:n) else as.vector(o)
  
  # Data frames
  results.efficiency <- matrix(NA, nrow = n, ncol = 1)
  results.iteration  <- matrix(NA, nrow = n, ncol = 1) 
  
  # Hyperplane consistency checking function
  on_h <- function(g){
    r.h         <- dm.sf(xdata, ydata, rts, g, wd = NULL, se, sg, date, cv, k)
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
    on.h        <- isTRUE(round(abs(dm.sf(z.xdata, z.ydata, rts, z.g, wd = NULL, T, sg, date, cv, k)$eff[k]), 9) == 0)
    new.g       <- abs(c(xdata[k,] - theta * xdata[k,], ydata[k,]/theta - ydata[k,]))
    results     <- list(test = on.h, theta = theta, new.g = new.g)
    return(results)
  }
  
  # Data modification for weak-disposability
  if(!is.null(wd)){
    x.o   <- unique(xdata)
    xdata <- rbind(xdata, x.o)
    ydata <- rbind(ydata, matrix(0, nrow(x.o), s, dimnames = list(NULL, names(ydata))))
    date  <- if(!is.null(date)) rbind(date, matrix(1, nrow(x.o), 1, dimnames = list(NULL, names(date))))
  }
  
  # Run sf
  g.o <- cbind(xdata, ydata)
  r   <- dm.sf(xdata, ydata, rts, g.o, wd = NULL, se, sg, date, cv, o)
  for(k in o){
    if(round(r$eff[k], 8) == 0){
      results.efficiency[k] <- 1
      results.iteration[k]  <- 1
    }else{
      h    <- 0
      flag <- T
      while(flag){
        z       <- on_h(g.o)
        flag    <- !z$test
        theta   <- z$theta
        g.o[k,] <- z$new.g
        h       <- h + 1
      }
      results.efficiency[k] <- theta
      results.iteration[k]  <- h
    }
  }
  temp.f <- dm.sf(xdata[1:n,, drop = F], ydata[1:n,, drop = F], rts, g.o, wd, se, sg, date[1:n,, drop = F], cv, o)
  
  # Get results
  results <- list(eff = results.efficiency, lambda = temp.f$lambda, mu = temp.f$mu,
                  xslack = temp.f$xslack, yslack = temp.f$yslack, iteration = results.iteration)
  return(results)
}
