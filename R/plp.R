plp <-
function (x){
  m       <- dim(x)[1]
  n       <- dim(x)[2]
  control <- lp.control(x)
  cat(paste("Model name: ", name.lp(x), "\n", sep = ""))
  ans     <- matrix(0, m + 1, n)
  for(j in 1:n) {
    col        <- get.column(x, j)
    col$column -> ans[1 + col$nzrow, j]
  }
  type  <- get.type(x);      kind <- get.kind(x)
  type[type == "integer"        ] <- "Int"
  type[type == "real"           ] <- "Real"
  kind[kind == "standard"       ] <- "Std"
  kind[kind == "semi-continuous"] <- "S-C"
  bounds                          <- get.bounds(x)
  upper <- bounds$upper;    lower <- bounds$lower
  ans   <- format(rbind(dimnames(x)[[2]], ans, kind, type, upper, lower), justify = "right")
  sense <- ifelse(control$sense == "minimize", "Minimize", "Maximize")
  lhs   <- get.constr.value(x, side = "lhs")
  rhs   <- get.constr.value(x, side = "rhs")
  r.nm  <- format(c("", sense, dimnames(x)[[1]], "Kind", "Type", "Upper", "Lower"))
  const <- format(c("", "", get.constr.type(x), "", "", "", ""), justify = "right")
  rhs   <- format(c("", "", as.character(rhs), "", "", "", ""), justify = "right")
  p.lhs <- any(!is.infinite(lhs[is.element(get.constr.type(x, as.char = FALSE), c(1, 2))]))
  lhs   <- format(c("", "", as.character(lhs), "", "", "", ""), justify = "right")
  if(p.lhs){
    ans <- cbind(r.nm, lhs, const, ans, const, rhs)
  }else{
    ans <- cbind(r.nm, ans, const, rhs)
  }
  ans   <- apply(ans, 1, paste, collapse = "  ")
  ans   <- paste(ans, collapse = "\n")
  m.nm  <- paste("Model name: ", name.lp(x), "\n", sep = "")
  ans   <- paste(m.nm, ans, "\n", sep = "")
  cat(ans)
  invisible(x)
}
