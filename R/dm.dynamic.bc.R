dm.dynamic.bc <-
function(xdata, ydata, zdata, bdata, rts = "crs", orientation = "i", wv = NULL){
  # Load library
  #library(lpSolveAPI)
  
  # Initial checks
  if(!identical(dim(xdata)[1], dim(ydata)[1], dim(zdata)[1]))                stop('Data must be balanced.')
  if(!identical(rev(dim(xdata))[1], rev(dim(ydata))[1], rev(dim(zdata))[1])) stop('Data must be balanced.')
  if(!(rts %in% c("crs", "vrs", "irs", "drs")))                              stop('rts must be "crs", "vrs", "irs", or "drs".')
  if(!(orientation %in% c("i", "o")))                                        stop('orientation must be "i", "o", or "n".')
  if(!all(wv >= 0))                                                          stop('wv must be >= 0.')
  
  # Parameters
  names <- if(is.null(rownames(xdata))) 1:n else rownames(xdata)
  xdata <- if(length(dim(xdata)) != 3) array(xdata, c(dim(xdata)[1], 1, dim(xdata)[2])) else as.array(xdata)
  ydata <- if(length(dim(ydata)) != 3) array(ydata, c(dim(ydata)[1], 1, dim(ydata)[2])) else as.array(ydata)
  zdata <- if(length(dim(zdata)) != 3) array(zdata, c(dim(zdata)[1], 1, dim(zdata)[2])) else as.array(zdata)
  bdata <- if(length(dim(bdata)) != 2) array(bdata, c(length(bdata), 1))                else as.array(bdata)
  n     <- dim(xdata)[1]
  m     <- dim(xdata)[2]
  s     <- dim(ydata)[2]
  b     <- dim(zdata)[2]
  t     <- dim(xdata)[3]
  wv    <- if(is.null(wv)) rep(1, t) else as.vector(wv)
  
  # Budget available at each T
  Z.A   <- array(bdata, c(n, b, t), dimnames = list(names, paste0("Avail.bg.", 1:b)))
  for(l in 1:b){for(k in 2:t){Z.A[, l, k] <- Z.A[, l, k - 1] - zdata[, l, (k - 1), drop = F]}}
  
  # Data frames
  results.lambda       <- array(NA, dim = c(n, n, t), dimnames = list(names, names))
  results.efficiency.s <- array(NA, dim = c(n, 1),    dimnames = list(names, "Eff.Sys"))
  results.efficiency.t <- array(NA, dim = c(n, t),    dimnames = list(names, paste0("Eff.T.", 1:t)))
  results.xslack       <- array(NA, dim = c(n, m, t), dimnames = list(names, paste0("xslack.", 1:m)))
  results.zslack       <- array(NA, dim = c(n, b, t), dimnames = list(names, paste0("zslack.", 1:b)))
  results.yslack       <- array(NA, dim = c(n, s, t), dimnames = list(names, paste0("yslack.", 1:s)))
  results.aslack       <- array(NA, dim = c(n, b, t), dimnames = list(names, paste0("aslack.", 1:b)))
  
  # Pointers
  p.eff <- n*t + 1
  p.xsl <- n*t + t + 1
  p.zsl <- n*t + t + m*t + 1
  p.ysl <- n*t + t + m*t + b*t + 1
  p.asl <- n*t + t + m*t + b*t + s*t + 1
  p.end <- n*t + t + m*t + b*t + s*t + b*t + 1
  
  # LP
  for(j in 1:n){
    
    # Declare LP
    lp.bc <- make.lp(0, n*t + t + m*t + b*t + s*t + b*t) # lambda + eff + xslack + zslack + yslack + aslack
    
    # Set objective
    if(orientation == "i") set.objfn(lp.bc,  wv, indices = c(p.eff:(p.xsl - 1)))
    if(orientation == "o") set.objfn(lp.bc, -wv, indices = c(p.eff:(p.xsl - 1)))
    
    # RTS
    if(rts == "vrs") for(k in 1:t){add.constraint (lp.bc, c(rep(1, n)), indices = c(((k-1)*n + 1):((k-1)*n + n)), "=", 1)}
    if(rts == "crs") for(k in 1:t){set.constr.type(lp.bc, 0, k)}
    if(rts == "irs") for(k in 1:t){add.constraint (lp.bc, c(rep(1, n)), indices = c(((k-1)*n + 1):((k-1)*n + n)), ">=", 1)}
    if(rts == "drs") for(k in 1:t){add.constraint (lp.bc, c(rep(1, n)), indices = c(((k-1)*n + 1):((k-1)*n + n)), "<=", 1)}
    
    # Each t
    for(k in 1:t){
      
      # Input constraint
      for(i in 1:m){
        if(orientation == "i"){
          add.constraint(lp.bc, c(xdata[, i, k], -xdata[j, i, k], 1), indices = c(((k-1)*n + 1):((k-1)*n + n), n*t + k, p.xsl - 1 + m*(k - 1) + i), "=", 0)
        }else{
          add.constraint(lp.bc, c(xdata[, i, k], 1), indices = c(((k-1)*n + 1):((k-1)*n + n), p.xsl - 1 + m*(k - 1) + i), "=", xdata[j, i, k])
        }
      }
      
      # Budget-spent constraint
      for(l in 1:b){
        if(orientation == "i"){
          add.constraint(lp.bc, c(zdata[, l, k], -zdata[j, l, k], 1), indices = c(((k-1)*n + 1):((k-1)*n + n), n*t + k, p.zsl - 1 + b*(k - 1) + l), "=", 0)
        }else{
          add.constraint(lp.bc, c(zdata[, l, k], 1), indices = c(((k-1)*n + 1):((k-1)*n + n), p.zsl - 1 + b*(k - 1) + l), "=", zdata[j, l, k])
        }
      }
      
      # Output constraint
      for(r in 1:s){
        if(orientation == "i"){
          add.constraint(lp.bc, c(ydata[, r, k], -1), indices = c(((k-1)*n + 1):((k-1)*n + n), p.ysl - 1 + s*(k - 1) + r), "=", ydata[j, r, k])
        }else{
          add.constraint(lp.bc, c(ydata[, r, k], -ydata[j, r, k], -1), indices = c(((k-1)*n + 1):((k-1)*n + n), n*t + k, p.ysl - 1 + s*(k - 1) + r), "=", 0)
        }
      }
      
      # Budget-available constraints
      for(l in 1:b){
        add.constraint(lp.bc, c(Z.A[, l, k], 1), indices = c(((k-1)*n + 1):((k-1)*n + n), p.asl - 1 + b*(k - 1) + l), "=", Z.A[j, l, k])  
      }
    }
    
    # Bounds
    set.bounds(lp.bc, lower = c(rep(0, n*t), rep(-Inf, t), rep(0, p.end - p.xsl)))
    
    # Solve
    solve.lpExtPtr(lp.bc)
    
    # Get results
    results.efficiency.s[j]  <- abs(get.objective(lp.bc))/sum(wv)
    temp.p                   <- get.variables(lp.bc)
    results.lambda[j,,]      <- array(temp.p[1:(n*t)], c(n, t))
    results.efficiency.t[j,] <- temp.p[p.eff:(p.xsl - 1)]
    results.xslack[j,,]      <- array(temp.p[p.xsl:(p.zsl - 1)], c(m, t))
    results.zslack[j,,]      <- array(temp.p[p.zsl:(p.ysl - 1)], c(b, t))
    results.yslack[j,,]      <- array(temp.p[p.ysl:(p.asl - 1)], c(s, t))
    results.aslack[j,,]      <- array(temp.p[p.asl:(p.end - 1)], c(1, t))
    
    # Stage II
    # Link previous solutions
    for(k in 1:t){add.constraint(lp.bc, c(1), indices = c(p.eff + k - 1), "=", results.efficiency.t[j, k])}

    # Slack sum max
    set.objfn(lp.bc, c(rep(-1, (p.end - p.xsl))), indices = c(p.xsl:(p.end - 1)))
    
    # solve
    solve.lpExtPtr(lp.bc)
    
    # Get results
    temp.s              <- get.variables(lp.bc)
    results.lambda[j,,] <- array(temp.s[1:(n*t)],           c(n, t))
    results.xslack[j,,] <- array(temp.s[p.xsl:(p.zsl - 1)], c(m, t))
    results.zslack[j,,] <- array(temp.s[p.zsl:(p.ysl - 1)], c(b, t))
    results.yslack[j,,] <- array(temp.s[p.ysl:(p.asl - 1)], c(s, t))
    results.aslack[j,,] <- array(temp.s[p.asl:(p.end - 1)], c(b, t))
  }
  
  # Store results
  results <- list(eff.s  = results.efficiency.s, 
                  eff.t  = results.efficiency.t, 
                  lambda = results.lambda, 
                  xslack = results.xslack, 
                  yslack = results.yslack, 
                  zslack = results.zslack,
                  aslack = results.aslack)
  return(results)
}
