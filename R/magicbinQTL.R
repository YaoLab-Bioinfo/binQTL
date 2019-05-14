

magicbinQTL <- function(phenotype="", genotype="") {
  kinship <- magicKinship(genotype)
  phe <- phenotype
  geno <- genotype
  
  mixed <- function(x,y,kk,method="REML",eigen=FALSE){
    loglike <- function(theta) {
      lambda <- exp(theta)
      logdt <- sum(log(lambda * delta + 1))
      h <- 1 / (lambda * delta + 1)
      yy <- sum(yu * h * yu)
      yx <- matrix(0, s, 1)
      xx <- matrix(0, s, s)
      for (i in 1:s) {
        yx[i] <- sum(yu * h * xu[, i])
        for (j in 1:s) {
          xx[i, j] <- sum(xu[, i] * h * xu[, j])
        }
      }
      xx
      if (method == "REML") {
        loglike <- -0.5 * logdt - 0.5 * (n - s) * log(yy - t(yx) %*% solve(xx) %*% yx) - 0.5 * log(det(xx))
      } else {
        loglike <- -0.5 * logdt - 0.5 * n * log(yy - t(yx) %*% solve(xx) %*% yx)
      }
      return(-loglike)
    }
    
    fixed <- function(lambda) {
      h <- 1 / (lambda * delta + 1)
      yy <- sum(yu * h * yu)
      yx <- matrix(0, s, 1)
      xx <- matrix(0, s, s)
      for (i in 1:s) {
        yx[i] <- sum(yu * h * xu[, i])
        for (j in 1:s) {
          xx[i, j] <- sum(xu[, i] * h * xu[, j])
        }
      }
      beta <- solve(xx, yx)
      if (method == "REML") {
        sigma2 <- (yy - t(yx) %*% solve(xx) %*% yx) / (n - s)
      } else {
        sigma2 <- (yy - t(yx) %*% solve(xx) %*% yx) / n
      }
      var <- diag(solve(xx) * drop(sigma2))
      stderr <- sqrt(var)
      return(c(beta, stderr, sigma2))
    }
    
    n <- length(y)
    qq <- eigen(kk, symmetric = TRUE)
    delta <- qq[[1]]
    uu <- qq[[2]]
    s <- ncol(x)
    yu <- t(uu) %*% y
    xu <- t(uu) %*% x
    theta <- 0
    parm <- optim(par = theta, fn = loglike, NULL, hessian = TRUE, method = "Brent", lower = -10, upper = 10)
    lambda <- exp(parm$par)
    conv <- parm$convergence
    fn1 <- parm$value
    fn0 <- loglike(-Inf)
    lrt <- 2 * (fn0 - fn1)
    hess <- parm$hessian
    parmfix <- fixed(lambda)
    beta <- parmfix[1:s]
    stderr <- parmfix[(s + 1):(2 * s)]
    ve <- parmfix[2 * s + 1]
    lod <- lrt / 4.61
    p_value <- 1 - pchisq(lrt, 1)
    va <- lambda * ve
    h2 <- va / (va + ve)
    par <- data.frame(method, beta, stderr, va, ve, lambda, h2, conv, fn1, fn0, lrt, lod, p_value)
    if(eigen){
      return(list(par,qq))
    } else {
      return(list(par))
    }
  }
  
  kk <- as.matrix(kinship)
  y <- as.matrix(phe[, 2])
  
  x <- matrix(1, length(y), 1)
  vy <- var(y)
  
  parm <- mixed(x, y, kk, eigen = T)
  lambda <- parm[[1]]$va[1] / parm[[1]]$ve[1]
  qq <- parm[[2]]
  delta <- qq[[1]]
  uu <- qq[[2]]
  
  yu <- t(uu) %*% y
  xu0 <- t(uu) %*% x
  xu00 <- cbind(xu0, y)
  
  fixed.parm <- function(lambda) {
    h <- 1 / (lambda * delta + 1)
    yy <- sum(yu * h * yu)
    yx <- matrix(0, w, 1)
    xx <- matrix(0, w, w)
    for (i in 1:w) {
      yx[i] <- sum(yu * h * xu[, i])
      for (j in 1:w) {
        xx[i, j] <- sum(xu[, i] * h * xu[, j])
      }
    } 
    beta <- solve(xx + diag(0.000000001, w), yx)
    sigma2 <- (yy - t(yx) %*% ginv(xx) %*% yx) / (n - w)
    var <- diag(ginv(xx) * drop(sigma2))
    stderr <- sqrt(var)
    gk <- beta[(w - s + 1):w]
    vk <- var[(w - s + 1):w]
    wk <- gk ^ 2 / (vk + 1e-8)
    beta <- beta[1:(w - s)]
    stderr <- stderr[1:(w - s)]
    value <- list(cbind(beta, stderr), c(gk, vk, wk), sigma2)
    return(value)
  }
  
  q <- ncol(xu00)
  
  name.fixed <- NULL
  name.error <- NULL
  for (i in 1:(q - 1)) {
    name.fixed <- c(name.fixed, paste("fixed_", i, sep = ""))
    name.error <- c(name.error, paste("error_", i, sep = ""))
  }
  
  geno <- as.matrix(genotype[, -c(1:4)])
  m <- nrow(geno)
  n <- ncol(geno)
  
  parr <- numeric()
  m1 <- m/8
  for(k in 1:m1) {
    zk <- t(geno[(8 * k - 7):(8 * k),])
    s <- ncol(zk)
    xu <- cbind(xu00[, -q], t(uu) %*% zk)
    w <- ncol(xu)
    parm <- fixed.parm(lambda)
    fixed <- t(parm[[1]][, 1])
    error <- t(parm[[1]][, 2])
    
    ve <- parm[[3]]
    va <- lambda * ve
    b1 <- parm[[2]][1:s]
    v1 <- parm[[2]][(s + 1):(2 * s)]
    w1 <- parm[[2]][(2 * s + 1):(3 * s)]
    wald1 <- sum(w1)
    p1 <- 1 - pchisq(wald1, (s - 1))
    lp1 <- -log10(p1)
    df <- max(s - 1 - sum(v1 / b1 ^ 2), 0) #change df
    pve <- t(b1) %*% var(zk) %*% b1 / vy
    t <- c(fixed, b1)
    par <- data.frame(pve, s, wald1, p1, lp1)[1, ]
    parr <- rbind(parr, par)
  }
  
  geno.map <- unique(genotype[, 1:4])
  binqtl.res <- data.frame(geno.map, parr, stringsAsFactors = FALSE)
  binqtl.res <- binqtl.res[, c(1:4, 8:9)]
  names(binqtl.res)[5:6] <- c("p", "-logp")
  return(binqtl.res)
}

