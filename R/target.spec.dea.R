target.spec.dea <-
function(xdata, ydata, date = NULL, t = NULL, dt = NULL, dmu, et = "c",
                            alpha = NULL, beta = NULL, wv = NULL, rts = "crs", sg = "ssm",
                            ftype = "d", ncv = NULL, env = NULL, cv = "convex", bound = TRUE, pin = TRUE){
  
  # Initial checks
  if(is.null(date) & sg != "ssm")                               stop('sg must be "ssm" when date is NULL.')
  if(is.null(date) | is.null(t) | is.null(dt))                  mtype <- "sidea" else mtype <- "tidea"
  if(mtype == "tidea") if(t <= min(date))                       stop('t is earlier than dataset.')
  if(mtype == "tidea") if(max(date) < t)                        stop('t is later than dataset.')
  if(dmu > nrow(as.matrix(xdata)) | dmu < 1)                    stop('dmu must indicate one of data points in the set.')
  if(!xor(is.null(alpha), is.null(beta)))                       stop('Either alpha or beta must be defined.')
  if(is.null(alpha) & !is.null(beta))                           orientation <- "i" else orientation <- "o"
  if(orientation == "i" & !(length(wv) %in% c(0, ncol(xdata)))) stop('wv must have the same length with input.')
  if(orientation == "o" & !(length(wv) %in% c(0, ncol(ydata)))) stop('wv must have the same length with output.')
  if(is.na(match(rts, c("crs", "vrs", "irs", "drs"))))          stop('rts must be "crs", "vrs", "irs", or "drs".')
  if(is.na(match(sg, c("ssm", "max", "min"))))                  stop('sg must be "ssm", "max", or "min".')
  if(is.na(match(ftype, c("d", "s"))))                          stop('ftype must be either "d" or "s".')
  if(is.na(match(cv, c("convex", "fdh"))))                      stop('cv must be "convex" or "fdh".')
  
  # Parameters 1
  xdata <- as.matrix(xdata)
  ydata <- as.matrix(ydata)
  date  <- if(mtype == "tidea") as.matrix(date)
  env   <- if(!is.null(env)) as.matrix(env)
  alpha <- if(orientation == "i") matrix(NA, nrow = 1, ncol = ncol(xdata)) else matrix(alpha, 1)
  beta  <- if(orientation == "o") matrix(NA, nrow = 1, ncol = ncol(ydata)) else matrix(beta,  1)
  wv.i  <- if(is.null(wv)){if(length(alpha) == 1) 1 else as.vector(1 - xdata[dmu,]/sum(xdata[dmu,]))}else{as.vector(wv)}
  wv.o  <- if(is.null(wv)){if(length(beta ) == 1) 1 else as.vector(1 - ydata[dmu,]/sum(ydata[dmu,]))}else{as.vector(wv)}
  
  # PPS
  if(mtype == "tidea"){
    # RoC
    roc_t <- roc.dea(xdata, ydata, date, t, rts, orientation, sg, ftype, ncv, env, cv)
    et    <- ifelse(et == "c", roc_t$eff_t[dmu,], et)
    
    # Reference set at t + dt
    id_soa        <- which(round(roc_t$eff_t, 8) == 1)
    d_f           <- date[id_soa,, drop = F]
    roc_temp      <- roc_t$roc_local[id_soa,]
    roc_t$roc_avg -> roc_temp[is.na(roc_temp)]
    x_f           <- if(orientation == "i") xdata[id_soa,, drop = F] * ((1 / roc_temp)^dt) else xdata[id_soa,, drop = F]
    y_f           <- if(orientation == "i") ydata[id_soa,, drop = F] else ydata[id_soa,, drop = F] * roc_temp^dt
    env_f         <- if(!is.null(env)) as.matrix(env[id_soa,])
    
    # Experiment set
    d_e   <- rbind(d_f, date[dmu,])
    x_e   <- if(orientation == "i") rbind(x_f, xdata[dmu,]) else rbind(x_f, alpha)
    y_e   <- if(orientation == "i") rbind(y_f, beta) else rbind(y_f, ydata[dmu,])
    env_e <- if(!is.null(env)) rbind(env_f, env[dmu,])
  }else{
    # Calc efficiency
    dea_t <- dm.dea(xdata, ydata, rts, orientation, ncv = ncv, env = env, cv = cv)
    et    <- ifelse(et == "c", dea_t$eff[dmu,], et)

    # Experiment set
    id_soa <- which(round(dea_t$eff, 8) == 1)
    x_e    <- if(orientation == "i") rbind(xdata[id_soa,, drop = F], xdata[dmu,]) else rbind(xdata[id_soa,, drop = F], alpha)
    y_e    <- if(orientation == "i") rbind(ydata[id_soa,, drop = F], beta) else rbind(ydata[id_soa,, drop = F], ydata[dmu,])
    env_e  <- if(!is.null(env)) rbind(as.matrix(env[id_soa,]), env[dmu,])
  }
  
  # Parameters 2
  n   <- nrow(x_e)
  m   <- ncol(x_e)
  s   <- ncol(y_e)
  rts <- ifelse(cv == "fdh", "vrs", rts)
  ncv <- if(is.null(ncv)) matrix(0, ncol = m + s) else as.matrix(ncv)

  # Feasibility check
  if(isTRUE(bound)){
    if(orientation == "i"){
      x_l <- rbind(x_e[-n,, drop = F], xdata[dmu,] * et)
      y_l <- rbind(y_e[-n,, drop = F], ydata[dmu,])
      fb  <- ydata[dmu,] * (dm.dea(x_l, y_l, rts, "o", ncv = ncv, env = env_e, cv = cv, o = n)$eff[n])
      if(sum(beta > fb) > 0) stop(paste0('Beta(', paste(beta, collapse = ", "),
                                         ') is greater than feasible bound(',
                                         paste(round(fb, 4), collapse = ", "), ').'))
    }else{
      x_l <- rbind(x_e[-n,, drop = F], xdata[dmu,])
      y_l <- rbind(y_e[-n,, drop = F], ydata[dmu,] * et)
      fb  <- xdata[dmu,] * (dm.dea(x_l, y_l, rts, "i", ncv = ncv, env = env_e, cv = cv, o = n)$eff[n])
      if(sum(alpha < fb) > 0) stop(paste0('Alpha(', paste(alpha, collapse = ", "),
                                          ') is smaller than feasible bound(',
                                          paste(round(fb, 4), collapse = ", "), ').'))
    } 
  }
  
  # Data frames
  lambda <- matrix(NA, nrow = 1, ncol = n)
  xslack <- matrix(NA, nrow = 1, ncol = m) 
  yslack <- matrix(NA, nrow = 1, ncol = s) 
  
  # Declare LP
  lp.idea <- if(orientation == "i") make.lp(0, n + m + m + s) else make.lp(0, n + s + m + s)

  # Set objective
  if(orientation == "i") set.objfn(lp.idea, c( wv.i), indices = c((n + 1):(n + m)))
  if(orientation == "o") set.objfn(lp.idea, c(-wv.o), indices = c((n + 1):(n + s)))
  
  # RTS
  if(rts == "vrs") add.constraint (lp.idea, c(rep(1, n)), indices = c(1:n), "=", 1)
  if(rts == "crs") set.constr.type(lp.idea, 0, 1)
  if(rts == "irs") add.constraint (lp.idea, c(rep(1, n)), indices = c(1:n), ">=", 1)
  if(rts == "drs") add.constraint (lp.idea, c(rep(1, n)), indices = c(1:n), "<=", 1)
  
  # Set type
  if(cv == "fdh") set.type(lp.idea, 1:n, "binary")
  
  # Input constraints
  for(i in 1:m){
    if(orientation == "i") add.constraint(lp.idea, c(x_e[, i], -et^(1 - ncv[1, i]), 1), indices = c(1:n, n + i, n + m + i), "=", 0)
    if(orientation == "o") add.constraint(lp.idea, c(x_e[, i], 1), indices = c(1:n, n + s + i), "=", alpha[1, i])
  }
  
  # Output constraints
  for(r in 1:s){
    if(orientation == "i") add.constraint(lp.idea, c(y_e[, r], -1), indices = c(1:n, n + m + m + r), "=", beta[1, r])
    if(orientation == "o") add.constraint(lp.idea, c(y_e[, r], -et^(1 - ncv[1, i]), -1), indices = c(1:n, n + r, n + s + m + r), "=", 0)
  }
  
  # External NDF
  if(!is.null(env_e)){
    for(j in 1:n){
      if(env_e[j, 1] < env_e[n, 1]){
        add.constraint(lp.idea, c(1), indices = c(j), "=", 0)
      }
    }
  }
  
  # Prevent perturbed DMU from being efficient
  if(!isTRUE(pin)){
    add.constraint(lp.idea, 1, indices = n, "=", 0)
  }
  
  # Bounds
  if(isTRUE(bound)){
    if(orientation == "i") set.bounds(lp.idea, upper = c(rep(Inf, n), xdata[dmu,], rep(Inf, m + s)))
    if(orientation == "o") set.bounds(lp.idea, lower = c(rep(0, n), ydata[dmu,], rep(0, m + s)))
  }

  # Solve
  solve.lpExtPtr(lp.idea)
  
  # Get results
  temp.p     <- get.variables(lp.idea)
  alpha [1,] <- if(orientation == "i") temp.p[(n + 1):(n + m)] else alpha
  beta  [1,] <- if(orientation == "o") temp.p[(n + 1):(n + s)] else beta

  # Stage II
  # Link previous solutions
  if(orientation == "i"){
    for(i in 1:m){
      add.constraint(lp.idea, c(1), indices = c(n + i), "=", alpha[1, i])
    }
  }else{
    for(r in 1:s){
      add.constraint(lp.idea, c(1), indices = c(n + r), "=", beta[1, r])
    }  
  }

  # Reset objective
  if(sg == "max") set.objfn(lp.idea, c(-d_e[1:n,]), indices = c(1:n))
  if(sg == "min") set.objfn(lp.idea, c( d_e[1:n,]), indices = c(1:n))
  if(sg == "ssm"){
    if(orientation == "i") set.objfn(lp.idea, c(rep(-1, m + s)), indices = c((n + m + 1):(n + m + m + s)))
    if(orientation == "o") set.objfn(lp.idea, c(rep(-1, m + s)), indices = c((n + s + 1):(n + s + m + s)))
  }
  
  # solve
  solve.lpExtPtr(lp.idea)
  
  # get results
  temp.s     <- get.variables(lp.idea)
  lambda[1,] <- temp.s[1:n]
  if(orientation == "i"){
    xslack[1,] <- temp.s[(n + m + 1):(n + m + m)]
    yslack[1,] <- temp.s[(n + m + m + 1):(n + m + m + s)]
  }else{
    xslack[1,] <- temp.s[(n + s + 1):(n + s + m)]
    yslack[1,] <- temp.s[(n + s + m + 1):(n + s + m + s)]  
  }

  results <- list(alpha = alpha, beta = beta, lambda = lambda, xslack = xslack, yslack = yslack)
  return(results)    
}
