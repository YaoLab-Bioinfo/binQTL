
genoConv <- function(genotype="") {
  ### Z matrix
  geno.z <- genotype[,-c(1:4)]
  geno.z[geno.z == "P1"] <- 1
  geno.z[geno.z == "P2"] <- -1
  geno.z[geno.z == "H"] <- 0
  geno.z <- as.matrix(geno.z)
  mode(geno.z) <- "numeric"

  z.mat <- t(geno.z) %*% geno.z
  geno.z <- cbind(genotype[, 1:4], geno.z)

  ### W matrix
  geno.w <- genotype[,-c(1:4)]
  geno.w[geno.w == "P1"] <- 0
  geno.w[geno.w == "P2"] <- 0
  geno.w[geno.w == "H"] <- 1
  geno.w <- as.matrix(geno.w)
  mode(geno.w) <- "numeric"

  w.mat <- t(geno.w) %*% geno.w
  geno.w <- cbind(genotype[, 1:4], geno.w)

  return(list(
    a.mat = z.mat,
    genotype.a = geno.z,
    d.mat = w.mat,
    genotype.d = geno.w
  ))
}
