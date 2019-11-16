target.arrival.sf <-
function(xdata, ydata, date, t, rts = "crs", g = NULL,
                              wd = NULL, sg = "ssm", ftype = "d", cv = "convex"){
  
  # Initial checks
  if(t <= min(date))                                   stop('t is earlier than dataset.')
  if(max(date) < t)                                    stop('t is later than dataset.')
  if(is.na(match(rts, c("crs", "vrs", "irs", "drs")))) stop('rts must be "crs", "vrs", "irs", or "drs".')
  if(is.na(match(sg,  c("ssm", "max", "min"))))        stop('sg must be "ssm", "max", or "min".')
  if(is.na(match(ftype, c("d","s"))))                  stop('ftype must be either "d" or "s".')
  if(is.na(match(cv,  c("convex", "fdh"))))            stop('cv must be "convex" or "fdh".')
  
  # Parameters
  xdata <- as.matrix(xdata)
  ydata <- as.matrix(ydata)
  date  <- as.matrix(date)
  g     <- if(is.null(g)) cbind(xdata, ydata) else as.matrix(g)
  n     <- nrow(xdata)
  m     <- ncol(xdata)
  s     <- ncol(ydata)
  wd    <- if(is.null(wd)) matrix(c(0), ncol = s) else as.matrix(wd)
  rts   <- ifelse(cv == "fdh", "vrs", rts)
  r     <- tail(which(sort(date) <= t), 1)
  id_f  <- which(date > t)
  
  # Data frames
  eff_t       <- array(NA, c(n, 1))
  eff_t_gm    <- array(NA, c(n, 1)) # Geometric mean for equi-ratio
  lambda      <- array(NA, c(n, n))
  ed          <- array(NA, c(n, 1))
  roc_ind     <- array(NA, c(n, 1))
  arrival_avg <- array(NA, c(n, 1))
  arrival_seg <- array(NA, c(n, 1))

  # Loop for eff_t
  for(i in id_f){
    # DMU set index
    id_s <- c(which(date <= t), i)
    
    # Run                        
    sf_t                 <- dm.sf(xdata[id_s,], ydata[id_s,], rts, g[id_s,], wd,
                                  1, sg, date[id_s,], cv, r + 1)
    eff_t[i,]            <- sf_t$eff[r + 1,]
    lambda[i, id_s[1:r]] <- sf_t$lambda[r + 1, 1:r]
  }
  
  # Effective date
  ed[id_f,] <- if(ftype == "s") rep(t, n - r) else lambda[id_f, -id_f] %*% date[-id_f,] / rowSums(lambda[id_f, -id_f])
  
  # RoC
  sf_roc         <- roc.sf(xdata, ydata, date, t, rts, g, wd, sg, ftype, cv)
  roc_local      <- sf_roc$roc_local
  roc_avg        <- sf_roc$roc_avg
  id_lroc        <- which(!is.na(roc_local))
  roc_ind[id_f,] <- lambda[id_f, id_lroc, drop = F] %*% roc_local[id_lroc,] / rowSums(lambda[id_f, id_lroc, drop = F])

  # Arrival target
  eff_t_gm[id_f,]      <- ((1 - eff_t[id_f,]) / (1 + eff_t[id_f,]))^0.5 # Geometric mean for equi-ratio
  roc_avg              -> roc_ind[!is.na(roc_ind) & roc_ind == 0 | is.nan(roc_ind),]  # replace 0 or NaN with roc_avg
  id_eff               <- id_f[is.finite(eff_t_gm[id_f]) & round(eff_t_gm[id_f], 8) > 1] 
  arrival_avg[id_eff,] <- ed[id_eff,] + log(eff_t_gm[id_eff,], exp(1)) / log(roc_avg, exp(1))
  arrival_seg[id_eff,] <- ed[id_eff,] + log(eff_t_gm[id_eff,], exp(1)) / log(roc_ind[id_eff,], exp(1))
  
  results <- list(eff_t = eff_t, lambda_t = lambda, eft_date = ed, 
                  roc_avg = roc_avg, roc_local = roc_local, roc_ind = roc_ind, 
                  arrival_avg = arrival_avg, arrival_seg = arrival_seg)
  return(results)    
}
