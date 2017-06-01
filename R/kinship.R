
kinship <- function(geno="") {
  z <- as.matrix(geno[, -c(1:4)])
  m <- nrow(z)
  n <- ncol(z)
  
  kk <- matrix(0, n, n)
  for (k in 1:m) {
    x <- as.factor(z[k, ])
    w <- model.matrix( ~ x - 1)
    kk <- kk + w %*% t(w)
  }
  kk <- kk / m
  kk <- data.frame(kk)
  
  colnames(kk) <- colnames(z)
  return(kk)
}
