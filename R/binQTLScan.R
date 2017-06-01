

binQTLScan <- function(phenotype="", genotype="", population="RIL") {
  Z <- as.matrix(genotype[, -c(1:4)])
  Y <- as.matrix(phenotype[, 2]) 
  genotype.n <- genoConv(genotype)
  kk <- kinship(genotype.n$genotype.a)
  m <- nrow(Z)
  n <- ncol(Z)
  X <- as.matrix(rep(1, n))

  qq <- eigen(kk, symmetric = T)
  uu <- qq$vectors
  delta <- qq$values

  r <- ncol(X)
  x <- t(uu) %*% X
  y <- t(uu) %*% Y
  if (population=="RIL") {
    L <- matrix(c(-0.5, 0.5), 1, 2)
  } else {
    L <- rbind(c(-0.5, 0, 0.5), c(-0.5, 1, -0.5))
  }
  invLL <- solve(L %*% t(L))
  rnk <- nrow(L)

  triplet <- TRUE
  parr <- numeric()
  for (k in 1:m) {
    if ((triplet == FALSE) | (k == 1) | (k == m)) {
      U <- as.factor(Z[k, ])
      Zk <- model.matrix( ~ U - 1)
      zk <- t(uu) %*% as.matrix(Zk)
      fit <- random( x = x, z = zk, y = y, m = m, kk = qq, cov = "qq", delta = delta, n = n)
    } else{
      U1 <- as.factor(Z[k - 1, ])
      U2 <- as.factor(Z[k, ])
      U3 <- as.factor(Z[k + 1, ])
      Z.1 <- model.matrix( ~ U1 - 1)
      Z.2 <- model.matrix( ~ U2 - 1)
      Z.3 <- model.matrix( ~ U3 - 1)
      Zk <- cbind(Z.1, Z.2, Z.3)
      zk <- t(uu) %*% as.matrix(Zk)
      fit <- random1( x = x, z = zk, y = y, m = m, kk = qq, cov = "qq", delta = delta, n = n )
    }
    covparm <- fit[[1]]
    vg <- covparm[1]
    vk <- covparm[2]
    ve <- covparm[3]
    lambda <- covparm[4]
    lambda.k <- covparm[5]
    beta <- fit[[2]]
    gamma <- fit[[3]]
    stderr <- fit[[4]]
    cov <- fit[[5]]
    lrt <- fit[[6]]
    effect <- L %*% gamma
    vv <- L %*% cov %*% t(L)
    wald <- t(effect) %*% solve(vv) %*% effect
    dk <- rnk - sum(diag(invLL %*% vv)) / vk
    if (dk < 1e-5) { pwald <- 1 } else { pwald <- 1 - pchisq(wald, dk) }
    if (lrt < 1e-8) { plrt <- 1 } else { plrt <- 0.5 * (1 - pchisq(lrt, 1)) }
    par <- c( k, vg, vk, ve, lambda, lambda.k, drop(beta), drop(gamma), drop(stderr),
        drop(effect), wald, dk, pwald, -log10(pwald), lrt, plrt, -log10(plrt)  )
    parr <- rbind(parr, par)
  }
  
  if (population=="RIL") {
    varnames <- c("k", "vg", "vk", "ve", "lambda", "lambdak", "beta", "g11",
      "g22", "sd11", "sd22")
    varnames <- c(varnames, "additive", "wald", "dk", "pwald", "logpwald",
      "lrt", "plrt", "logplrt")
  } else {
    varnames <- c("k", "vg", "vk", "ve", "lambda", "lambdak", "beta", "g11",
                  "g12", "g22", "sd11", "sd12", "sd22")
    varnames <- c(varnames, "additive", "dominance", "wald", "dk", "pwald",
                  "logpwald", "lrt", "plrt", "logplrt")
  }
  
  out <- as.matrix(parr)
  colnames(out) <- varnames
  rownames(out) <- 1:nrow(out)
  result <- data.frame(genotype[, 1:4], out, stringsAsFactors=FALSE)
  result <- result[, c("Bin", "Chr", "Start", "Stop", "plrt", "logplrt")]
  names(result)[5:6] <- c("p", "-logp")
  return(result)
}
