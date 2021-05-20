dm.network.dea <-
function(xdata.s1, ydata.s1 = NULL, zdata, xdata.s2 = NULL, ydata.s2, 
                           rts = "crs", orientation = "i", type = "nc", leader = "1st", ss = 10^-4, o = NULL){
  
  # Initial checks
  if(is.na(match(rts, c("crs", "vrs", "irs", "drs")))) stop('rts must be "crs", "vrs", "irs", or "drs".')
  if(is.na(match(orientation, c("i", "o"))))           stop('orientation must be either "i" or "o".')
  if(is.na(match(type, c("co", "nc"))))                stop('type must be "co" or "nc".')
  if(is.na(match(leader, c("1st", "2nd"))))            stop('leader must be "1st" or "2nd".')
  if(!is.null(o) && !all(o <= nrow(xdata.s1)))         stop('o must be element(s) of n.')
  
  # Parameters
  xdata.s1 <- as.matrix(xdata.s1)
  ydata.s2 <- as.matrix(ydata.s2)
  zdata    <- as.matrix(zdata)
  if(!is.null(ydata.s1)){ydata.s1 <- as.matrix(ydata.s1)}
  if(!is.null(xdata.s2)){xdata.s2 <- as.matrix(xdata.s2)}  
  n        <- nrow(xdata.s1)
  m.s1     <- ncol(xdata.s1)
  s.s1     <- ifelse(is.null(ydata.s1), 0, ncol(ydata.s1))
  p        <- ncol(zdata)
  m.s2     <- ifelse(is.null(xdata.s2), 0, ncol(xdata.s2))
  s.s2     <- ncol(ydata.s2)
  o        <- if(is.null(o)) c(1:n) else as.vector(o)
  
  # Data frames
  res.eff.s1 <- array(NA, c(n, 1))
  res.eff.s2 <- array(NA, c(n, 1))
  res.v.s1   <- array(NA, c(n, m.s1)) 
  res.u.s1   <- array(NA, c(n, s.s1)) 
  res.p      <- array(NA, c(n, p)) 
  res.w.s1   <- array(NA, c(n, 1)) 
  res.v.s2   <- array(NA, c(n, m.s2)) 
  res.u.s2   <- array(NA, c(n, s.s2)) 
  res.w.s2   <- array(NA, c(n, 1))

  # Indices for coding convenience
  no.dv.t <- (m.s1 + s.s1 + p + 1 + m.s2 + s.s2 + 1) 
  id.v.s1 <- 1:m.s1
  id.u.s1 <- if(is.null(ydata.s1)) 0 else id.v.s1 + m.s1
  id.p    <- (m.s1 + s.s1 + 1):(m.s1 + s.s1 + p)
  id.w.s1 <-  m.s1 + s.s1 + p + 1
  id.v.s2 <- if(is.null(xdata.s2)) 0 else (m.s1 + s.s1 + p + 2):(m.s1 + s.s1 + p + 1 + m.s2)
  id.u.s2 <- (m.s1 + s.s1 + p + 1 + m.s2 + 1):(m.s1 + s.s1 + p + 1 + m.s2 + s.s2)
  id.w.s2 <- no.dv.t
  
  # Analysis
  if(type == "co"){ # Centralized / Cooperative game
    if(!is.null(ydata.s1) | !is.null(xdata.s2)){ # With external input/output
      for(k in o){
        # Declare LP
        lp.ndea <- make.lp(0, no.dv.t) # v1+u1+p+w1+v2+u2+w2
        
        # Labeling
        temp <- c(paste0("v1", 1:m.s1), if(!is.null(ydata.s1)) paste0("u1", 1:s.s1),
                  paste0("p", 1:p), paste0("w1"), if(!is.null(xdata.s2)) paste0("v2", 1:m.s2), 
                  paste0("u2", 1:s.s2), paste0("w2"))
        dimnames(lp.ndea)[[2]] <- temp
        
        # Objective
        if(orientation == "o") stop('DJL is too busy to implement OO model.')
        if(orientation == "i") set.objfn(lp.ndea,  c(if(is.null(ydata.s1)) NULL else -ydata.s1[k,], -zdata[k,], 1), 
                                         indices = c(if(is.null(ydata.s1)) NULL else id.u.s1, id.p, id.w.s1))
        
        # RTS
        if(rts == "crs") add.constraint(lp.ndea, c(1, 1), indices = c(id.w.s1, id.w.s2), "=", 0)
        if(rts == "irs") stop('DJL is too busy to implement IRS model.')
        if(rts == "drs") stop('DJL is too busy to implement DRS model.')
        
        # Constraint for o
        add.constraint(lp.ndea, xdata.s1[k,], indices = id.v.s1, "=", 1)
        
        # Constraint for all
        for(d in o){
          # Stage 1
          add.constraint(lp.ndea, c(-xdata.s1[d,], 
                                    if(is.null(ydata.s1)) NULL else ydata.s1[d,],
                                    zdata[d,], -1), 
                         indices = c(id.v.s1, 
                                     if(is.null(ydata.s1)) NULL else id.u.s1, 
                                     id.p, id.w.s1), "<=", 0)
          
          # Stage 2
          add.constraint(lp.ndea, c(-zdata[d,],
                                    if(is.null(xdata.s2)) NULL else -xdata.s2[d,],
                                    ydata.s2[d,], -1), 
                         indices = c(id.p, 
                                     if(is.null(xdata.s2)) NULL else id.v.s2, 
                                     id.u.s2, id.w.s2), "<=", 0)  
        }
        
        # Bounds
        temp.lb <- rep(0, no.dv.t)
        if(rts == "vrs"){temp.lb[c(id.w.s1, id.w.s2)] <- -Inf}
        set.bounds(lp.ndea, lower = temp.lb)  
        
        # Solve
        solve.lpExtPtr(lp.ndea)
        
        # Get maximum possible results of Stage 1
        res.eff.max <- abs(get.objective(lp.ndea))
        
        # Heuristic search
        res.eff.cand <- c()
        for(q in 1:floor(res.eff.max/ss)){
          # Running status 
          if(q == ceiling(res.eff.max/ss*0.25)){
            print(paste("Heuristic search:  25% of job done for DMU", k))  
          }else if(q == ceiling(res.eff.max/ss*0.5)){
            print(paste("Heuristic search:  50% of job done for DMU", k))  
          }else if(q == ceiling(res.eff.max/ss*0.75)){
            print(paste("Heuristic search:  75% of job done for DMU", k))  
          }else if(q == floor(res.eff.max/ss)){
            print(paste("Heuristic search: 100% of job done for DMU", k))  
          }

          # Temporal eff
          res.eff.temp <- res.eff.max - (q - 1) * ss
          
          # Declare LP
          lp.ndea <- make.lp(0, no.dv.t) # v1+u1+p+w1+v2+u2+w2
          
          # Labeling
          temp <- c(paste0("v1", 1:m.s1), if(!is.null(ydata.s1)) paste0("u1", 1:s.s1),
                    paste0("p", 1:p), paste0("w1"), if(!is.null(xdata.s2)) paste0("v2", 1:m.s2), 
                    paste0("u2", 1:s.s2), paste0("w2"))
          dimnames(lp.ndea)[[2]] <- temp
          
          # Objective
          if(orientation == "o") stop('DJL is too busy to implement OO model.')
          if(orientation == "i") set.objfn(lp.ndea,  c(-ydata.s2[k,] * res.eff.temp, 1 * res.eff.temp), 
                                           indices = c(id.u.s2, id.w.s2))
          
          # RTS
          if(rts == "crs") add.constraint(lp.ndea, c(1, 1), indices = c(id.w.s1, id.w.s2), "=", 0)
          if(rts == "irs") stop('DJL is too busy to implement IRS model.')
          if(rts == "drs") stop('DJL is too busy to implement DRS model.')
          
          # Constraint for o
          add.constraint(lp.ndea, c(zdata[k,], if(is.null(xdata.s2)) NULL else xdata.s2[k,]), 
                         indices = c(id.p, if(is.null(xdata.s2)) NULL else id.v.s2), "=", 1)
          
          # Eff linking constraint
          add.constraint(lp.ndea, c(-xdata.s1[k,] * res.eff.temp, 
                                    if(is.null(ydata.s1)) NULL else ydata.s1[k,], zdata[k,], -1), 
                         indices = c(id.v.s1, if(is.null(ydata.s1)) NULL else id.u.s1, id.p, id.w.s1), "=", 0)
          
          # Constraint for all
          for(d in o){
            # Stage 1
            add.constraint(lp.ndea, c(-xdata.s1[d,], 
                                      if(is.null(ydata.s1)) NULL else ydata.s1[d,],
                                      zdata[d,], -1), 
                           indices = c(id.v.s1, 
                                       if(is.null(ydata.s1)) NULL else id.u.s1, 
                                       id.p, id.w.s1), "<=", 0)
            
            # Stage 2
            add.constraint(lp.ndea, c(-zdata[d,],
                                      if(is.null(xdata.s2)) NULL else -xdata.s2[d,],
                                      ydata.s2[d,], -1), 
                           indices = c(id.p, 
                                       if(is.null(xdata.s2)) NULL else id.v.s2, 
                                       id.u.s2, id.w.s2), "<=", 0)  
          }
          
          # Bounds
          temp.lb <- rep(0, no.dv.t)
          if(rts == "vrs"){temp.lb[c(id.w.s1, id.w.s2)] <- -Inf}
          set.bounds(lp.ndea, lower = temp.lb)  
          
          # Solve
          solve.lpExtPtr(lp.ndea)
          
          # Get results
          res.eff.cand <- c(res.eff.cand, abs(get.objective(lp.ndea)))
        }
        
        # Pick the best system eff
        id.best <- which(res.eff.cand == max(res.eff.cand))
        
        # Obtain eff.s2 with the best system eff
        # Fixate eff.s1
        res.eff.s1[k,] <- res.eff.max - (id.best - 1) * ss
        
        # Declare LP
        lp.ndea <- make.lp(0, no.dv.t) # v1+u1+p+w1+v2+u2+w2
        
        # Labeling
        temp <- c(paste0("v1", 1:m.s1), if(!is.null(ydata.s1)) paste0("u1", 1:s.s1),
                  paste0("p", 1:p), paste0("w1"), if(!is.null(xdata.s2)) paste0("v2", 1:m.s2), 
                  paste0("u2", 1:s.s2), paste0("w2"))
        dimnames(lp.ndea)[[2]] <- temp
        
        # Objective
        if(orientation == "o") stop('DJL is too busy to implement OO model.')
        if(orientation == "i") set.objfn(lp.ndea,  c(-ydata.s2[k,] * res.eff.s1[k,], 1 * res.eff.s1[k,]), 
                                         indices = c(id.u.s2, id.w.s2))
        
        # RTS
        if(rts == "crs") add.constraint(lp.ndea, c(1, 1), indices = c(id.w.s1, id.w.s2), "=", 0)
        if(rts == "irs") stop('DJL is too busy to implement IRS model.')
        if(rts == "drs") stop('DJL is too busy to implement DRS model.')
        
        # Constraint for o
        add.constraint(lp.ndea, c(zdata[k,], if(is.null(xdata.s2)) NULL else xdata.s2[k,]), 
                       indices = c(id.p, if(is.null(xdata.s2)) NULL else id.v.s2), "=", 1)
        
        # Eff linking constraint
        add.constraint(lp.ndea, c(-xdata.s1[k,] * res.eff.s1[k,], 
                                  if(is.null(ydata.s1)) NULL else ydata.s1[k,], zdata[k,], -1), 
                       indices = c(id.v.s1, if(is.null(ydata.s1)) NULL else id.u.s1, id.p, id.w.s1), "=", 0)
        
        # Constraint for all
        for(d in o){
          # Stage 1
          add.constraint(lp.ndea, c(-xdata.s1[d,], 
                                    if(is.null(ydata.s1)) NULL else ydata.s1[d,],
                                    zdata[d,], -1), 
                         indices = c(id.v.s1, 
                                     if(is.null(ydata.s1)) NULL else id.u.s1, 
                                     id.p, id.w.s1), "<=", 0)
          
          # Stage 2
          add.constraint(lp.ndea, c(-zdata[d,],
                                    if(is.null(xdata.s2)) NULL else -xdata.s2[d,],
                                    ydata.s2[d,], -1), 
                         indices = c(id.p, 
                                     if(is.null(xdata.s2)) NULL else id.v.s2, 
                                     id.u.s2, id.w.s2), "<=", 0)  
        }
        
        # Bounds
        temp.lb <- rep(0, no.dv.t)
        if(rts == "vrs"){temp.lb[c(id.w.s1, id.w.s2)] <- -Inf}
        set.bounds(lp.ndea, lower = temp.lb)  
        
        # Solve
        solve.lpExtPtr(lp.ndea)
        
        # Get results
        res.all         <- get.variables(lp.ndea)
        res.v.s1[k,]    <- res.all[id.v.s1]
        res.u.s1[k,]    <- if(is.null(ydata.s1)) NA else res.all[id.u.s1]
        res.p[k,]       <- res.all[id.p]
        res.w.s1[k,]    <- res.all[id.w.s1]
        res.v.s2[k,]    <- if(is.null(xdata.s2)) NA else res.all[id.v.s2]
        res.u.s2[k,]    <- res.all[id.u.s2]
        res.w.s2[k,]    <- res.all[id.w.s2]
        res.eff.s2.temp <- sum(res.u.s2[k,] * ydata.s2[k,]) - res.w.s2[k,] - sum(c(res.p[k,], res.v.s2[k,]) * c(zdata[k,], xdata.s2[k,]))
        res.eff.s2[k,]  <- ifelse(res.eff.s2.temp == 0, 1,
                                  (sum(res.u.s2[k,] * ydata.s2[k,]) - res.w.s2[k,]) / (sum(c(res.p[k,], res.v.s2[k,]) * c(zdata[k,], xdata.s2[k,]))))

      }

    }else{ # No external input/output
      for(k in o){
        # Declare LP
        lp.ndea <- make.lp(0, no.dv.t) # v1+u1+p+w1+v2+u2+w2
        
        # Labeling
        temp <- c(paste0("v1", 1:m.s1), if(!is.null(ydata.s1)) paste0("u1", 1:s.s1),
                  paste0("p", 1:p), paste0("w1"), if(!is.null(xdata.s2)) paste0("v2", 1:m.s2), 
                  paste0("u2", 1:s.s2), paste0("w2"))
        dimnames(lp.ndea)[[2]] <- temp
        
        # Objective
        if(orientation == "o") stop('DJL is too busy to implement OO model.')
        if(orientation == "i") set.objfn(lp.ndea, c(1, -ydata.s2[k,], 1), indices = c(id.w.s1, id.u.s2, id.w.s2))
        
        # RTS
        if(rts == "crs") add.constraint(lp.ndea, c(1, 1), indices = c(id.w.s1, id.w.s2), "=", 0)
        if(rts == "irs") stop('DJL is too busy to implement IRS model.')
        if(rts == "drs") stop('DJL is too busy to implement DRS model.')
        
        # Constraint for o
        add.constraint(lp.ndea, xdata.s1[k,], indices = id.v.s1, "=", 1)
        
        # Constraint for all
        for(d in o){
          # Stage 1
          add.constraint(lp.ndea, c(-xdata.s1[d,], 
                                    if(is.null(ydata.s1)) NULL else ydata.s1[d,],
                                    zdata[d,], -1), 
                         indices = c(id.v.s1, 
                                     if(is.null(ydata.s1)) NULL else id.u.s1, 
                                     id.p, id.w.s1), "<=", 0)
          
          # Stage 2
          add.constraint(lp.ndea, c(-zdata[d,],
                                    if(is.null(xdata.s2)) NULL else -xdata.s2[d,],
                                    ydata.s2[d,], -1), 
                         indices = c(id.p, 
                                     if(is.null(xdata.s2)) NULL else id.v.s2, 
                                     id.u.s2, id.w.s2), "<=", 0)  
        }
        
        # Bounds
        temp.lb <- rep(0, no.dv.t)
        if(rts == "vrs"){temp.lb[c(id.w.s1, id.w.s2)] <- -Inf}
        set.bounds(lp.ndea, lower = temp.lb)  
        
        # Solve
        solve.lpExtPtr(lp.ndea)
        
        # Get results
        res.all         <- get.variables(lp.ndea)
        res.v.s1[k,]    <- res.all[id.v.s1]
        res.u.s1[k,]    <- if(is.null(ydata.s1)) NA else res.all[id.u.s1]
        res.p[k,]       <- res.all[id.p]
        res.w.s1[k,]    <- res.all[id.w.s1]
        res.v.s2[k,]    <- if(is.null(xdata.s2)) NA else res.all[id.v.s2]
        res.u.s2[k,]    <- res.all[id.u.s2]
        res.w.s2[k,]    <- res.all[id.w.s2]
        res.eff.s1[k,]  <- sum(c(res.u.s1[k,], res.p[k,]) * c(ydata.s1[k,], zdata[k,])) - res.w.s1[k,]
        res.eff.s2.temp <- sum(res.u.s2[k,] * ydata.s2[k,]) - res.w.s2[k,] - sum(c(res.p[k,], res.v.s2[k,]) * c(zdata[k,], xdata.s2[k,]))
        res.eff.s2[k,]  <- ifelse(res.eff.s2.temp == 0, 1,
                                  (sum(res.u.s2[k,] * ydata.s2[k,]) - res.w.s2[k,]) / (sum(c(res.p[k,], res.v.s2[k,]) * c(zdata[k,], xdata.s2[k,]))))
      }
    }

  }else{ # Decentralized / Stackelberg game

    if(leader == "1st"){
      # Leader
      res.eff.s1 <- dm.dea(xdata.s1, cbind(ydata.s1, zdata), rts, orientation)$eff
      
      # Follower
      for(k in o){
        # Declare LP
        lp.ndea <- make.lp(0, no.dv.t) # v1+u1+p+w1+v2+u2+w2
        
        # Labeling
        temp <- c(paste0("v1", 1:m.s1), if(!is.null(ydata.s1)) paste0("u1", 1:s.s1),
                  paste0("p", 1:p), paste0("w1"), if(!is.null(xdata.s2)) paste0("v2", 1:m.s2), 
                  paste0("u2", 1:s.s2), paste0("w2"))
        dimnames(lp.ndea)[[2]] <- temp
        
        # Objective
        if(orientation == "o") stop('DJL is too busy to implement OO model.')
        if(orientation == "i") set.objfn(lp.ndea, c(-ydata.s2[k,], 1), indices = c(id.u.s2, id.w.s2))
        
        # RTS
        if(rts == "crs") add.constraint(lp.ndea, c(1, 1), indices = c(id.w.s1, id.w.s2), "=", 0)
        if(rts == "irs") stop('DJL is too busy to implement IRS model.')
        if(rts == "drs") stop('DJL is too busy to implement DRS model.')
        
        # Constraint for o
        add.constraint(lp.ndea, c(-xdata.s1[k,] * res.eff.s1[k,], ydata.s1[k,], zdata[k,], -1), 
                       indices = c(id.v.s1, if(is.null(ydata.s1)) NULL else id.u.s1, id.p, id.w.s1), "=", 0)
        
        add.constraint(lp.ndea, c(if(is.null(xdata.s2)) NULL else xdata.s2[k,], zdata[k,]), 
                       indices = c(if(is.null(xdata.s2)) NULL else id.v.s2, id.p), "=", 1)
        
        # Constraint for all
        for(d in o){
          # Stage 1
          add.constraint(lp.ndea, c(-xdata.s1[d,], 
                                    if(is.null(ydata.s1)) NULL else ydata.s1[d,],
                                    zdata[d,], -1), 
                         indices = c(id.v.s1, 
                                     if(is.null(ydata.s1)) NULL else id.u.s1, 
                                     id.p, id.w.s1), "<=", 0)
          
          # Stage 2
          add.constraint(lp.ndea, c(-zdata[d,],
                                    if(is.null(xdata.s2)) NULL else -xdata.s2[d,],
                                    ydata.s2[d,], -1), 
                         indices = c(id.p, 
                                     if(is.null(xdata.s2)) NULL else id.v.s2, 
                                     id.u.s2, id.w.s2), "<=", 0)  
        }
        
        # Bounds
        temp.lb <- rep(0, no.dv.t)
        if(rts == "vrs"){temp.lb[c(id.w.s1, id.w.s2)] <- -Inf}
        set.bounds(lp.ndea, lower = temp.lb)  
        
        # Solve
        solve.lpExtPtr(lp.ndea)
        
        # Plan B for multiple optima
        if(solve.lpExtPtr(lp.ndea) == 3){
          # Constraint for o: v1x1 = 1
          add.constraint(lp.ndea, xdata.s1[k,], indices = id.v.s1, "=", 1)

          # Re-solve
          solve.lpExtPtr(lp.ndea)
        }
          
        # Get results
        res.all        <- get.variables(lp.ndea)
        res.v.s1[k,]   <- res.all[id.v.s1]
        res.u.s1[k,]   <- if(is.null(ydata.s1)) NA else res.all[id.u.s1]
        res.p[k,]      <- res.all[id.p]
        res.w.s1[k,]   <- res.all[id.w.s1]
        res.v.s2[k,]   <- if(is.null(xdata.s2)) NA else res.all[id.v.s2]
        res.u.s2[k,]   <- res.all[id.u.s2]
        res.w.s2[k,]   <- res.all[id.w.s2]
        res.eff.s2[k,] <- sum(res.u.s2[k,] * ydata.s2[k,]) - res.w.s2[k,]
      }
      
    }else{
      # Leader
      res.eff.s2 <- dm.dea(cbind(xdata.s2, zdata), ydata.s2, rts, orientation)$eff
      
      # Follower
      for(k in o){
        # Declare LP
        lp.ndea <- make.lp(0, no.dv.t) # v1+u1+p+w1+v2+u2+w2
        
        # Labeling
        temp <- c(paste0("v1", 1:m.s1), if(!is.null(ydata.s1)) paste0("u1", 1:s.s1),
                  paste0("p", 1:p), paste0("w1"), if(!is.null(xdata.s2)) paste0("v2", 1:m.s2), 
                  paste0("u2", 1:s.s2), paste0("w2"))
        dimnames(lp.ndea)[[2]] <- temp
        
        # Objective
        if(orientation == "o") stop('DJL is too busy to implement OO model.')
        if(orientation == "i") set.objfn(lp.ndea, c(-zdata[k,], 
                                                    if(is.null(ydata.s1)) NULL else -ydata.s1[k,], 1),
                                         indices = c(id.p, 
                                                     if(is.null(ydata.s1)) NULL else id.u.s1, id.w.s1))
        
        # RTS
        if(rts == "crs") add.constraint(lp.ndea, c(1, 1), indices = c(id.w.s1, id.w.s2), "=", 0)
        if(rts == "irs") stop('DJL is too busy to implement IRS model.')
        if(rts == "drs") stop('DJL is too busy to implement DRS model.')
        
        # Constraint for o
        add.constraint(lp.ndea, c(-zdata[k,] * res.eff.s2[k,], 
                                  if(is.null(xdata.s2)) NULL else -xdata.s2[k,] * res.eff.s2[k,], 
                                  ydata.s2[k,], -1), 
                       indices = c(id.p, 
                                   if(is.null(xdata.s2)) NULL else id.v.s2, 
                                   id.u.s2, id.w.s1), "=", 0)

        add.constraint(lp.ndea, c(xdata.s1[k,]), indices = c(id.v.s1), "=", 1)
        
        # Constraint for all
        for(d in o){
          # Stage 1
          add.constraint(lp.ndea, c(-xdata.s1[d,], 
                                    if(is.null(ydata.s1)) NULL else ydata.s1[d,],
                                    zdata[d,], -1), 
                         indices = c(id.v.s1, 
                                     if(is.null(ydata.s1)) NULL else id.u.s1, 
                                     id.p, id.w.s1), "<=", 0)
          
          # Stage 2
          add.constraint(lp.ndea, c(-zdata[d,],
                                    if(is.null(xdata.s2)) NULL else -xdata.s2[d,],
                                    ydata.s2[d,], -1), 
                         indices = c(id.p, 
                                     if(is.null(xdata.s2)) NULL else id.v.s2, 
                                     id.u.s2, id.w.s2), "<=", 0)  
        }
        
        # Bounds
        temp.lb <- rep(0, no.dv.t)
        if(rts == "vrs"){temp.lb[c(id.w.s1, id.w.s2)] <- -Inf}
        set.bounds(lp.ndea, lower = temp.lb)  
        
        # Solve
        solve.lpExtPtr(lp.ndea)
        
        # Plan B for multiple optima
        if(solve.lpExtPtr(lp.ndea) == 3){
          # Constraint for o: v1x1 = 1
          add.constraint(lp.ndea, xdata.s1[k,], indices = id.v.s1, "=", 1)
          
          # Re-solve
          solve.lpExtPtr(lp.ndea)
        }
        
        # Get results
        res.all        <- get.variables(lp.ndea)
        res.v.s1[k,]   <- res.all[id.v.s1]
        res.u.s1[k,]   <- if(is.null(ydata.s1)) NA else res.all[id.u.s1]
        res.p[k,]      <- res.all[id.p]
        res.w.s1[k,]   <- res.all[id.w.s1]
        res.v.s2[k,]   <- if(is.null(xdata.s2)) NA else res.all[id.v.s2]
        res.u.s2[k,]   <- res.all[id.u.s2]
        res.w.s2[k,]   <- res.all[id.w.s2]
        res.eff.s1[k,] <- sum(c(res.u.s1[k,], res.p[k,]) * c(ydata.s1[k,], zdata[k,])) - res.w.s1[k,]
      }
    }
    
  }
  
  # Returning object
  results <- list(eff.s1 = res.eff.s1, eff.s2 = res.eff.s2, 
                  v.s1 = res.v.s1, u.s1 = res.u.s1, p = res.p, w.s1 = res.w.s1,
                  v.s2 = res.v.s2, u.s2 = res.u.s2, w.s2 = res.w.s2)
  
  return(results)
}
