
aovQTL <- function(phenotype="", genotype="") {
  genotype.n <- genoConv(genotype)
  genotype.a <- genotype.n$genotype.a
  genotype.d <- genotype.n$genotype.d
  z <- t(genotype.a[, -c(1:4)])
  w <- t(genotype.d[, -c(1:4)])
  y <- as.matrix(phenotype[, 2])
  m <- ncol(z)

  #Simple least square method
  parr.lst <- lapply(k<-1:m, function(k) {
    zk <- cbind(z[, k], w[, k])
    p <- anova(lm(y ~ zk))[1, 5]

    return(c(k, p))
  })

  parr.df <- do.call(rbind, parr.lst)
  res.return <- data.frame(genotype.a[, 1:4], parr.df[, 2])
  res.return$LOD <- -log10(res.return[,5])
  names(res.return) <- c(names(genotype.a)[1:4], "p", "-logp")

  return(res.return)
}

