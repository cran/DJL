roc.malmquist <-
function(xdata, ydata, tm = NULL, dm = "dea", rts = "crs", orientation = "n", g = NULL, 
                          wd = NULL, ncv = NULL, env = NULL, cv = "convex"){
  
  # Initial checks
  if(!(3 %in% c(length(dim(xdata)), length(dim(ydata)))))              stop('Data must be 3-dimensional.')
  if(dim(xdata)[length(dim(xdata))] != dim(ydata)[length(dim(ydata))]) stop('Data must be balanced.')
  if(!is.null(tm) && length(tm) != dim(xdata)[length(dim(xdata))])     stop('tm must have a length of the time horizon.')
  if(dm == "dea" && orientation == "n")                                stop('DEA orientation must be either "i", or "o".')
  if(is.na(match(dm,          c("dea", "ddf", "hdf", "sbm", "sf"))))   stop('dm must be "dea", "ddf", "hdf", "sbm", or "sf".')
  if(is.na(match(rts,         c("crs", "vrs", "irs", "drs"))))         stop('rts must be "crs", "vrs", "irs", or "drs".')
  if(is.na(match(orientation, c("i", "o", "n"))))                      stop('orientation must be "i", "o", or "n".')
  if(is.na(match(cv,          c("convex", "fdh"))))                    stop('cv must be "convex" or "fdh".')
  
  # Parameters
  xdata <- if(length(dim(xdata)) != 3) array(xdata, c(dim(xdata)[1], 1, dim(xdata)[2])) else as.array(xdata)
  ydata <- if(length(dim(ydata)) != 3) array(ydata, c(dim(ydata)[1], 1, dim(ydata)[2])) else as.array(ydata)
  g.exe <- paste0("array(c(", 
                  paste0("xdata[,,", 1:dim(xdata)[3], "],", "ydata[,,", 1:dim(ydata)[3], "]", collapse=","), 
                  "), c(dim(xdata)[1], dim(xdata)[2] + dim(ydata)[2], dim(xdata)[3]))")
  g     <- if(is.null(g)) eval(parse(text = g.exe)) else as.array(g)
  n     <- dim(xdata)[1]
  m     <- dim(xdata)[2]
  s     <- dim(ydata)[2]
  wd    <- if(is.null(wd)) matrix(c(0), ncol = s) else as.matrix(wd)
  rts   <- ifelse(cv == "fdh", "vrs", rts)
  t     <- ifelse(is.null(tm), dim(xdata)[length(dim(xdata))], length(tm))
  tm    <- if(is.null(tm)) paste0("t", 1:t) else as.vector(tm)
  
  # Model Arguments
  if(dm %in% c("dea", "sbm")){
    m.arg <- list(rts = rts, orientation = orientation)
  }else if(dm == "hdf"){
    m.arg <- list(rts = rts, wd = wd)
  }else{
    m.arg <- list(rts = rts, g = g, wd = wd)
  }
  
  # Inter-temporal dm
  inter <- function(f, t){
    temp.it <- vector()
    for(j in 1:n){
      temp.x  <- rbind(xdata[j,, f], as.matrix(xdata[,, t]))
      temp.y  <- rbind(ydata[j,, f], as.matrix(ydata[,, t]))
      temp.g  <- if(is.null(g)) cbind(temp.x, temp.y) else 
        temp.se <- do.call(paste0("dm.", dm), append(list(xdata = temp.x, ydata = temp.y, se = T, o = 1), m.arg))
      temp.no <- do.call(paste0("dm.", dm), append(list(xdata = temp.x, ydata = temp.y, se = F, o = 1), m.arg))
      temp.it <- c(temp.it, ifelse(round(temp.no$eff[1], 5) < 1, temp.no$eff[1], temp.se$eff[1]))
    }
    return(temp.it)
  }
  
  # Malmquist Index
  cu <- fs <- mi <- data.frame()
  for(i in 1:(t-1)){
    # Within efficiency
    m.0.0 <- do.call(paste0("dm.", dm), append(list(xdata = xdata[,, i],     ydata = ydata[,, i]),     m.arg))$eff
    m.1.1 <- do.call(paste0("dm.", dm), append(list(xdata = xdata[,, i + 1], ydata = ydata[,, i + 1]), m.arg))$eff
    
    # Inter-temporal efficiency
    m.0.1 <- inter(i,     i + 1)
    m.1.0 <- inter(i + 1, i    )
    
    # Tick mark
    temp.tm <- rep(paste0(tm[i], "-", tm[i + 1]), n)    
    
    # Catching-Up (CU aka TEC) Index
    temp.cu <- m.1.1 / m.0.0
    cu      <- rbind(cu, data.frame(Period = temp.tm, DMU = factor(1:n), CU = temp.cu))
    
    # Frontier-Shift (FS) Index
    temp.fs <- ((m.0.0 / m.0.1) * (m.1.0 / m.1.1))^0.5
    fs      <- rbind(fs, data.frame(Period = temp.tm, DMU = factor(1:n), FS = temp.fs))
    
    # Malmquist Index (MI)
    temp.mi <- temp.cu * temp.fs
    mi      <- rbind(mi, data.frame(Period = temp.tm, DMU = factor(1:n), MI = temp.mi))
  }
  
  # Retuning results object
  results <- list(cu = cu, fs = fs, mi = mi)
  return(results)
}
