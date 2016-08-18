map.corr <-
function(data, from = "median", threshold = 0.3, r.name = FALSE){
  
  # Initial checks
  if(is.na(match(from, c("mean", "median")))) stop('"from" must be either "mean" or "median".')
  if(dim(data)[2] != 2)                       stop('"data" must be a bivariate data.frame/matrix.')

  # Parameter
  k     <- 0
  n     <- nrow(data)
  r.eff <- sum(complete.cases(data))
  
  # Run
  results <- matrix(NA, n + 1, 3, dimnames = list(0 : n, c("corr", "drop", "n")))
  for(i in 0 : n){
    if(i == 0){
      data.eff <- data
      r.kill   <- NA
    }else{
      if(sum(!is.na(data.eff)) == 4 | all(round(cor(data.eff,use = "complete.obs"), 5) == 1)) break
      r.kill             <- dm.mahalanobis(data.eff, from = from, p = 100)$suspect[1]
      data.eff[r.kill, ] <- NA
    }
    
    if(!all(apply(data.eff[complete.cases(data.eff), ], 2, sd) > 0)) break
    results[i + 1, 1] <- cor(data.eff[, 1], data.eff[, 2], use = "complete.obs")
    results[i + 1, 2] <- r.kill
    results[i + 1, 3] <- ifelse(i == 0, r.eff, sum(complete.cases(data.eff)))
    if(i == 0){
      next
    }else if(abs(results[i, 1] - results[i + 1, 1]) > threshold){
      k <- k + 1
    }
  }
  
  # Detect odd jump
  if(k > 0){
    jump <- array(NA, c(2, 2, k))
    q    <- 0
    for(p in 1 : (i - 1)){
      if(abs(results[p, 1] - results[p + 1, 1]) > threshold){
        jump[, , q + 1] <- results[p : (p + 1), c(1, 3)]
        q               <- q + 1
        if(q == k) break
      }    
    }
  }
  
  # Parameter
  gradient.b <- colorRampPalette(c("red", "lightsalmon", "whitesmoke", "cornflowerblue", "blue"))
  gradient.d <- colorRampPalette(c("maroon", "darksalmon", "darkgray", "steelblue", "darkblue"))
  c.mean     <- round(mean(results[, 1], na.rm = T), 2) * 100 + 101
  c.total    <- round(results[1, 1], 2) * 100 + 101
  c.point    <- round(results[1 : i, 1], 2) * 100 + 101
  m          <- mean(results[, 1], na.rm = T)
  m.s        <- ifelse(m < 0, "-", ifelse(m > 0, "+", ""))
  t          <- results[1, 1]
  t.s        <- ifelse(t < 0, "-", ifelse(t > 0, "+", ""))
  cex.c      <- ifelse(n < 50, 2, 1.3)
  cex.t      <- ifelse(n < 50, 1, 0.7)

  # Plot
  plot(results[, 3], results[, 1], 
       ylim = c(-1.3, 1.3), 
       xlab = "Data point (n)", ylab = "Correlation (r)", 
       main = "Correlation Reliability Test", yaxt = "n")
  
  axis(2, -1,   "-1.0", col.axis = "red",            las = 2)
  axis(2, -0.5, "-0.5", col.axis = "lightsalmon",    las = 2)
  axis(2, 0,    "0.0",  col.axis = "darkgray",       las = 2)
  axis(2, 0.5,  "0.5",  col.axis = "cornflowerblue", las = 2)
  axis(2, 1,    "1.0",  col.axis = "blue",           las = 2)
  
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "whitesmoke", lty = 0)
  grid(NULL, NULL, lty = 1, lwd = 2, col = "white") 
  lines(results[order(results[, 3]), 3], results[order(results[, 3]), 1], 
        col = "dimgray", lwd = 0.5)
  
  if(k > 0){
    for(d in 1 : k){
      lines(jump[, 2, d], jump[, 1, d], col = "lawngreen", lwd = 2)
      text(mean(jump[, 2, d]), mean(jump[, 1, d]) + 0.2, "!", 
           col = "lawngreen", cex = cex.t + 1)
    }
  }
  
  points(results[, 3], results[, 1], 
         pch = 21, cex = cex.c, lwd = 0.3, 
         col = "dimgray", bg = gradient.b(201)[c.point])
  text(mean(results[, 3], na.rm = T),  1.2, 
       paste0("Total Correlation: ",   t.s,sprintf("%.2f", abs(round(t, 2)))), 
       col = gradient.d(201)[c.total])
  text(mean(results[, 3], na.rm = T), -1.2, 
       paste0("Average Correlation: ", m.s,sprintf("%.2f", abs(round(m, 2)))), 
       col = gradient.d(201)[c.mean])
  text(results[2 : i, 3], results[2 : i, 1] - 0.1, 
       if(r.name) rownames(data)[results[2 : i, 2]] else results[2 : i, 2], 
       cex = cex.t, 
       col = "dimgray")
  
  # Return results
  if(r.name) results[1 : i, 2] <- c("None", rownames(data)[results[2 : i, 2]]) else results[1, 2] <- "None"
  reliability <- data.frame(results[1 : i, ])
  results     <- list(reliability = reliability)
  return(results) 
}
