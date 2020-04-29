roc.dea <-
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
  o     <- matrix(c(1:n), ncol = 1) # original data order
  r     <- tail(which(sort(date) <= t), 1)
  
  # Sort data ascending order
  x   <- xdata[order(date),, drop = F]
  y   <- ydata[order(date),, drop = F]
  d   <- date [order(date),, drop = F]
  o   <- o    [order(date),, drop = F]
  env <- if(!is.null(env)) env[order(date),, drop = F]
  
  # Data frames
  eff_r     <- array(NA, c(n, 1))
  eff_t     <- array(NA, c(n, 1))
  lambda    <- array(NA, c(n, n))
  ed        <- array(NA, c(n, 1))
  sl        <- array(NA, c(n, 1))
  roc       <- array(NA, c(n, 1))
  local_roc <- array(NA, c(n, 1))
  
  # Loop for eff_r and eff_t
  for(i in unique(d[1:r])){
    # Run                        
    dea_r                 <- dm.dea(subset(x, d <= i), subset(y, d <= i), rts, orientation,
                                    0, sg, subset(d, d <= i), ncv, env, cv, which(d == i))
    eff_r[which(d == i),] <- dea_r$eff[which(d == i),]
    if(i == d[r]){
      dea_t               <- dm.dea(subset(x, d <= i), subset(y, d <= i), rts, orientation,
                                    0, sg, subset(d, d <= i), ncv, env, cv)
      eff_t[1:r,]         <- dea_t$eff[1:r,]
      lambda[1:r, 1:r]    <- dea_t$lambda[1:r, 1:r]
    } 
  }
  
  # Effective date
  ed <- if(ftype == "s") rep(t, r) else lambda[, 1:r] %*% d[1:r,] / rowSums(lambda, na.rm=T)
  
  # RoC
  id_roc       <- which(round(eff_r[, 1], 8) == 1 & round(eff_t[, 1], 8) != 1 & ed[, 1] > d[, 1])
  delta_t      <- 1/(ed - d)
  roc[id_roc,] <- if(orientation == "i") (1 / eff_t[id_roc,])^delta_t[id_roc,] else (eff_t[id_roc,])^delta_t[id_roc,]
  
  # RoC filter
  roc[roc[id_roc,] > 10,] <- NA
  avgroc                  <- mean(roc, na.rm = T)
  
  # RoC segmentation
  id_local_roc             <- which(colSums(lambda, na.rm = T) > 0)
  temp                     <- t(lambda[id_roc, id_local_roc, drop = F]) %*% roc[id_roc] / colSums(lambda[id_roc, id_local_roc, drop = F])
  temp[is.nan(temp),]      <- NA # For coding convenience, could be improved
  local_roc[id_local_roc,] <- temp

  # Sort results back to original order
  eff_r     <- eff_r[order(o),,           drop = F]
  eff_t     <- eff_t[order(o),,           drop = F]
  lambda    <- lambda[order(o), order(o), drop = F]
  ed        <- ed[order(o),,              drop = F]
  roc       <- roc[order(o),,             drop = F]
  local_roc <- local_roc[order(o),,       drop = F]
  
  results <- list(eff_r = eff_r, eff_t = eff_t, lambda_t = lambda, eft_date = ed,
                  roc_past = roc, roc_local = local_roc, roc_avg = avgroc)
  return(results)
}
