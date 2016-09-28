target.arrival.dea <-
function(xdata, ydata, date, t, rts = "crs", orientation,
                               sg = "ssm", ftype = "d", ncv = NULL, env = NULL, cv = "convex"){
  
  # Initial checks
  if(is.na(match(rts, c("crs", "vrs", "irs", "drs")))) stop('rts must be "crs", "vrs", "irs", or "drs".')
  if(is.na(match(orientation, c("i", "o"))))           stop('orientation must be either "i" or "o".')
  if(is.na(match(sg, c("ssm", "max", "min"))))         stop('sg must be "ssm", "max", or "min".')
  if(is.na(match(ftype,c("d","s"))))                   stop('ftype must be either "d" or "s".')
  if(t <= min(date))                                   stop('t is earlier than dataset.')
  if(max(date) < t)                                    stop('t is later than dataset.')
  if(is.na(match(cv, c("convex", "fdh"))))             stop('cv must be "convex" or "fdh".')
  
  # Parameters
  xdata <- as.matrix(xdata)
  ydata <- as.matrix(ydata)
  date  <- as.matrix(date)
  env   <- if(!is.null(env)) as.matrix(env)
  n     <- nrow(xdata)
  m     <- ncol(xdata)
  s     <- ncol(ydata)
  rts   <- ifelse(cv == "fdh", "vrs", rts)
  ncv   <- if(is.null(ncv)) matrix(0, ncol = m + s) else as.matrix(ncv)
  env   <- if(!is.null(env)) as.matrix(env)
  r     <- tail(which(sort(date) <= t), 1)
  id_f  <- which(date > t)
  
  # Data frames
  eff_t       <- array(NA, c(n, 1))
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
    dea_t                <- dm.dea(xdata[id_s,], ydata[id_s,], rts, orientation,
                                   1, sg, date[id_s,], ncv, env[id_s,], cv, r + 1)
    eff_t[i,]            <- dea_t$eff[r + 1,]
    lambda[i, id_s[1:r]] <- dea_t$lambda[r + 1, 1:r]
  }
  
  # Effective date
  ed[id_f,] <- if(ftype == "s") rep(t, n - r) else lambda[id_f, -id_f] %*% date[-id_f,] / rowSums(lambda[id_f, -id_f])
  
  # RoC
  roc_t          <- roc.dea(xdata, ydata, date, t, rts, orientation, sg, ftype, ncv, env, cv)
  roc_local      <- roc_t$roc_local
  roc_avg        <- roc_t$roc_avg
  id_lroc        <- which(!is.na(roc_local))
  roc_ind[id_f,] <- lambda[id_f, id_lroc] %*% roc_local[id_lroc,] / rowSums(lambda[id_f, id_lroc])
    
  # Arrival target
  eff_temp             <- if(orientation == "i") eff_t else 1/eff_t # for coding convenience
  roc_avg              -> roc_ind[!is.na(roc_ind) & roc_ind == 0 | is.nan(roc_ind),]  # replace 0 or NaN with roc_avg
  id_eff               <- id_f[is.finite(eff_t[id_f]) & round(eff_temp[id_f], 8) > 1] 
  arrival_avg[id_eff,] <- ed[id_eff,] + log(eff_temp[id_eff,], exp(1)) / log(roc_avg, exp(1))
  arrival_seg[id_eff,] <- ed[id_eff,] + log(eff_temp[id_eff,], exp(1)) / log(roc_ind[id_eff,], exp(1))
  
  results <- list(eff_t = eff_t, lambda_t = lambda, eft_date = ed, roc_avg = roc_avg, roc_local = roc_local,
                  roc_ind = roc_ind, arrival_avg = arrival_avg, arrival_seg = arrival_seg)
  return(results)    
}
