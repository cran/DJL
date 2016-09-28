dm.sbm <-
function(xdata, ydata, rts = "crs", 
                   orientation = "n", se = FALSE, sg = "ssm", date = NULL, cv = "convex", o = NULL){

  # Initial checks
  if(is.na(match(rts, c("crs", "vrs", "irs", "drs")))) stop('rts must be "crs", "vrs", "irs", or "drs".')
  if(is.na(match(orientation, c("n", "i", "o"))))      stop('orientation must be "n", "i" or "o".')
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
  se    <- ifelse(is.logical(se), ifelse(isTRUE(se), 1, 0), se)
  rts   <- ifelse(cv == "fdh", "vrs", rts)
  o     <- if(is.null(o)) c(1:n) else as.vector(o)
  
  # Data frames
  results.efficiency <- matrix(NA, nrow = n, ncol = 1)
  results.lambda     <- matrix(NA, nrow = n, ncol = n)
  results.xslack     <- matrix(NA, nrow = n, ncol = m) 
  results.yslack     <- matrix(NA, nrow = n, ncol = s) 
  results.xtarget    <- matrix(NA, nrow = n, ncol = m) 
  results.ytarget    <- matrix(NA, nrow = n, ncol = s) 
  
  # LP
  if(se == 0){
    for(k in o){   
      # Declare LP
      lp.sbm <- make.lp(0, (n + 1 + m + s)) # lambda+t+xslack+yslack
      
      # Set objective
      if(orientation == "o") set.objfn(lp.sbm, c(-1, -1 / (s * ydata[k, ])),
                                       indices = c(n + 1, (n + 1 + m + 1):(n + 1 + m + s)))
      if(orientation != "o") set.objfn(lp.sbm,c(1, -1 / (m * xdata[k, ])), 
                                       indices=c((n + 1):(n + 1 + m)))
      
      # RTS
      if(rts == "vrs") add.constraint(lp.sbm, c(rep(1, n)), indices = c(1:n), "=", 1)
      if(rts == "crs") set.constr.type(lp.sbm, 0, 1)
      if(rts == "irs") add.constraint(lp.sbm, c(rep(1, n)), indices = c(1:n), ">=", 1)
      if(rts == "drs") add.constraint(lp.sbm, c(rep(1, n)), indices = c(1:n), "<=", 1)
      
      # Set type
      if(cv == "fdh") set.type(lp.sbm, 1:n, "binary")
      
      # Let t be 1 to avoid clutters for non-CRS / oriented model
      if(rts != "crs")       add.constraint(lp.sbm, c(1), indices = c(n + 1), "=", 1)
      if(orientation != "n") add.constraint(lp.sbm, c(1), indices = c(n + 1), "=", 1)
      
      # Normalization constraint
      if(orientation == "n") add.constraint(lp.sbm, c(1, 1 / (s * ydata[k,])),
                                            indices = c(n + 1, (n + 1 + m + 1):(n + 1 + m + s)), "=", 1)
      
      # Input constraints
      for(i in 1:m){
        add.constraint(lp.sbm, c(xdata[, i], -xdata[k, i], 1),
                       indices = c(1:(n + 1), n + 1 + i), "=", 0)
      }
      
      # Output constraints
      for(r in 1:s){
        add.constraint(lp.sbm, c(ydata[, r], -ydata[k, r], -1),
                       indices = c(1:(n + 1), n + 1 + m + r), "=", 0)
      }
      
      # Oriented model - no slacks
      if(orientation == "i") add.constraint(lp.sbm, rep(1, s), 
                                            indices = c((n + 1 + m + 1):(n + 1 + m + s)), "=", 0)
      if(orientation == "o") add.constraint(lp.sbm, rep(1, m),
                                            indices = c((n + 1 + 1):(n + 1 + m)), "=", 0)
      
      # Bounds
      set.bounds(lp.sbm, lower = c(rep(0, n + 1 + m + s)))  

      # Solve
      solve.lpExtPtr(lp.sbm)
      
      # Get results
      results.efficiency[k] <- ifelse(orientation == "o", -1 / get.objective(lp.sbm), get.objective(lp.sbm))
      temp.p                <- get.variables(lp.sbm)
      results.lambda[k,]    <- temp.p[1:n] / temp.p[n + 1]
      results.xslack[k,]    <- temp.p[(n + 2):(n + 1 + m)] / temp.p[n + 1]
      results.yslack[k,]    <- temp.p[(n + 1 + m + 1):(n + 1 + m + s)] / temp.p[n + 1]
      for(i in 1:m){results.xtarget[k, i] <- sum(results.lambda[k,] * xdata[, i])}
      for(r in 1:s){results.ytarget[k, r] <- sum(results.lambda[k,] * ydata[, r])}
      
      # Stage II
      if(exists("sg")){
        # Link previous solutions
        if(orientation == "o") add.constraint(lp.sbm, c(-1, -1 / (s * ydata[k,])),
                                              indices = c(n + 1, (n + 1 + m + 1):(n + 1 + m + s)),
                                              "=", results.efficiency[k])
        if(orientation != "o") add.constraint(lp.sbm, c(1, -1 / (m * xdata[k, ])), 
                                              indices = c((n + 1):(n + 1 + m)), 
                                              "=", results.efficiency[k])
        
        # date sum
        if(sg == "max") set.objfn(lp.sbm, c(-date[1:n]), indices = c(1:n))
        if(sg == "min") set.objfn(lp.sbm, c(date[1:n]),  indices = c(1:n))
        
        # slack sum max
        if(sg == "ssm") set.objfn(lp.sbm, c(rep(-1, m + s)), indices = c((n + 2):(n + 1 + m + s)))
        
        # solve
        solve.lpExtPtr(lp.sbm)
        
        # get results
        temp.s             <- get.variables(lp.sbm)
        results.lambda[k,] <- temp.s[1:n] / temp.s[n + 1]
        results.xslack[k,] <- temp.s[(n + 2):(n + 1 + m )] / temp.s[n + 1]
        results.yslack[k,] <- temp.s[(n + 1 + m + 1):(n + 1 + m + s)] / temp.s[n + 1]
        for(i in 1:m){results.xtarget[k, i] <- sum(results.lambda[k,] * xdata[, i])}
        for(r in 1:s){results.ytarget[k, r] <- sum(results.lambda[k,] * ydata[, r])}
      }
    }
  }
  if(se == 1){
    for(k in o){   
      # Declare LP
      lp.sbm <- make.lp(0, (n + 1 + m + s)) # lambda+t+xb+yb / xtarget<-xb, ytarget<-yb
      
      # Set objective
      if(orientation == "o") set.objfn(lp.sbm, c(-1 / (s * ydata[k,])),
                                       indices = c((n + 1 + m + 1):(n + 1 + m + s)))
      if(orientation != "o") set.objfn(lp.sbm, c(1 / (m * xdata[k,])),
                                       indices = c((n + 2):(n + 1 + m)))
      
      # RTS
      if(rts == "vrs") add.constraint(lp.sbm, c(rep(1, n)), indices = c(1:n), "=", 1)
      if(rts == "crs") set.constr.type(lp.sbm, 0, 1)
      if(rts == "irs") add.constraint(lp.sbm, c(rep(1, n)), indices = c(1:n), ">=", 1)
      if(rts == "drs") add.constraint(lp.sbm, c(rep(1, n)), indices = c(1:n), "<=", 1)
      
      # Set type
      if(cv == "fdh") set.type(lp.sbm, 1:n, "binary")
      
      # Let t be 1 to avoid clutters for non-CRS / oriented model
      if(rts != "crs")       add.constraint(lp.sbm, c(1), indices = c(n + 1), "=", 1)
      if(orientation != "n") add.constraint(lp.sbm, c(1), indices = c(n + 1), "=", 1)

      # Normalization constraint
      if(orientation == "n") add.constraint(lp.sbm, c(1 / (s * ydata[k,])),
                                            indices = c((n + 1 + m + 1):(n + 1 + m + s)), "=", 1)
      
      # Input constraints
      for(i in 1:m){
        add.constraint(lp.sbm, c(xdata[, i], -1),  indices = c(1:n, n + 1 + i), "<=", 0)
        add.constraint(lp.sbm, c(xdata[k, i], -1), indices = c(n + 1, n + 1 + i), "<=", 0)
        if(orientation == "o") add.constraint(lp.sbm, c(xdata[k, i], -1), 
                                              indices = c((n + 1), (n + 1 + i)), "=", 0)
      }
      
      # Output constraints
      for(r in 1:s){
        add.constraint(lp.sbm, c(ydata[, r], -1), indices = c(1:n, n + 1 + m + r), ">=", 0)
        add.constraint(lp.sbm, c(ydata[k,r],-1),  indices = c(n + 1, n + 1 + m + r), ">=", 0)
        if(orientation == "i") add.constraint(lp.sbm, c(ydata[k, r], -1),
                                              indices = c((n + 1), (n + 1 + m + r)), "=", 0)
      }
      
      # Reference set for super
      add.constraint(lp.sbm, c(1), indices = c(k), "=", 0)
      
      # Bounds
      set.bounds(lp.sbm, lower = c(rep(0, n + 1 + m + s)))

      # Solve
      solve.lpExtPtr(lp.sbm)
      
      # Get results
      results.efficiency[k] <- ifelse(orientation == "o", -1 / get.objective(lp.sbm), get.objective(lp.sbm))
      temp.p                <- get.variables(lp.sbm)
      results.lambda[k,]    <- temp.p[1:n] / temp.p[n + 1]
      results.xtarget[k,]   <- temp.p[(n + 2):(n + 1 + m)] / temp.p[n + 1]
      results.ytarget[k,]   <- temp.p[(n + 1 + m + 1):(n + 1 + m + s)] / temp.p[n + 1]
      for(i in 1:m){results.xslack[k, i] <- results.xtarget[k,i] - sum(results.lambda[k,] * xdata[, i])}
      for(r in 1:s){results.yslack[k, r] <- sum(results.lambda[k,] * ydata[, r]) - results.ytarget[k, r]}
      
      # Stage II
      if(exists("sg")){
        # Link previous solutions
        if(orientation == "o") add.constraint(lp.sbm, c(-1 / (s * ydata[k,])),
                                              indices = c((n + 1 + m + 1):(n + 1 + m + s)),
                                              "=", results.efficiency[k])
        if(orientation != "o") add.constraint(lp.sbm, c(1 / (m * xdata[k,])),
                                              indices = c((n + 2):(n + 1 + m)),
                                              "=", results.efficiency[k])
        
        # date sum
        if(sg == "max") set.objfn(lp.sbm, c(-date[1:n]), indices = c(1:n))
        if(sg == "min") set.objfn(lp.sbm, c(date[1:n]),  indices = c(1:n))
        
        # slack sum max
        if(sg == "ssm") set.objfn(lp.sbm, c(rep(-1, m + s)), indices = c((n + 2):(n + 1 + m + s)))
        
        # solve
        solve.lpExtPtr(lp.sbm)
        
        # get results
        temp.s              <- get.variables(lp.sbm)
        results.lambda[k,]  <- temp.s[1:n] / temp.s[n + 1]
        results.xtarget[k,] <- temp.s[(n + 2):(n + 1 + m)] / temp.s[n + 1]
        results.ytarget[k,] <- temp.s[(n + 1 + m + 1):(n + 1 + m + s)] / temp.s[n + 1]
        for(i in 1:m){results.xslack[k, i] <- results.xtarget[k, i] - sum(results.lambda[k,] * xdata[, i])}
        for(r in 1:s){results.yslack[k, r] <- sum(results.lambda[k,] * ydata[, r]) - results.ytarget[k, r]}
      }
    }
  }
  results <- list(eff = results.efficiency, lambda = results.lambda, xslack = results.xslack,
                  yslack = results.yslack, xtarget = results.xtarget, ytarget = results.ytarget)
  return(results)
}
