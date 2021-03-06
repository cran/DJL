map.soa.sbm <-
function(xdata, ydata, date, rts = "crs", orientation = "n",
                        sg = "ssm", cv = "convex", mk = "dmu"){
  
  # Initial checks
  if(is.na(match(rts, c("crs", "vrs", "irs", "drs")))) stop('rts must be "crs", "vrs", "irs", or "drs".')
  if(is.na(match(orientation, c("n", "i", "o"))))      stop('orientation must be "n", "i", or "o".')
  if(is.na(match(sg, c("ssm", "max", "min"))))         stop('sg must be "ssm", "max", or "min".')
  if(is.na(match(mk, c("dmu", "eff"))))                stop('mk must be either "dmu" or "eff".')
  if(is.na(match(cv, c("convex", "fdh"))))             stop('cv must be "convex" or "fdh".')
  
  # Parameters
  xdata <- as.matrix(xdata)
  ydata <- as.matrix(ydata)
  date  <- if(!is.null(date)) as.matrix(date)
  n     <- nrow(xdata)
  m     <- ncol(xdata)
  s     <- ncol(ydata)
  rts   <- ifelse(cv == "fdh", "vrs", rts)
  o     <- matrix(c(1:n), ncol = 1) # original data order
  ud    <- sort(unique(date))
  l     <- length(ud)
  
  # Sort data ascending order
  x <- xdata[order(date),, drop = F]
  y <- ydata[order(date),, drop = F]
  d <- date [order(date),, drop = F]
  o <- o    [order(date),, drop = F]
  
  # Map frame
  map.soa <- matrix(NA, n, l, dimnames = list(NULL, ud)) 
  
  # Generate the map
  for(i in ud){
    # run
    sbm.t <- dm.sbm(subset(x, d <= i), subset(y, d <= i), rts, orientation, 0,
                    sg, subset(d, d <= i), cv)
    
    # SOA index
    #id.soa <- which(round(sbm.t$eff, 8) == 1) # if slacks are not concerned
    id.soa <- which(round(sbm.t$eff, 8) == 1 & 
                    rowSums(cbind(round(sbm.t$xslack, 8), 
                                  round(sbm.t$yslack, 8))) == 0)
      
    # Mapping
    if(mk == "dmu"){
      if(i == ud[1]){
        map.soa[1:length(id.soa), 1] <- o[id.soa]
      }else{
        p <- which(ud == i)
        for(k in 1:length(id.soa)){
          id.preb <- which(map.soa[, p - 1] == o[id.soa[k],])
          if(length(id.preb) > 0){
            map.soa[id.preb, p] <- o[id.soa[k],]
          }else{
            map.soa[sum(rowSums(map.soa, na.rm = T) > 0) + 1, p] <- o[id.soa[k],]
          }
        }
      }
    }else{
      gsoa <- if(i == ud[1]) id.soa else union(gsoa, id.soa)
      map.soa[1:length(gsoa), which(ud == i)] <- sbm.t$eff[gsoa,]
    }
  }
  
  # Prune the map
  map.soa           <- map.soa[1:max(which(!is.na(map.soa[, l]))),] 
  rownames(map.soa) <- if(mk == "dmu") unique(na.omit(c(map.soa))) else c(o[gsoa,]) 
  
  # Print
  print(map.soa)
}
