

random <- function(x, z, y, m, kk = qq, cov = "qq", delta, n) {
  if (cov == "kk") {
    qq <- eigen(kk, symmetric = T)
    uu <- qq$vectors; delta <- qq$values
    x <- t(uu) %*% x; y <- t(uu) %*% y; z <- t(uu) %*% z
  }
  r <- ncol(x); q <- ncol(z)
  
  func <- function(par) {
    lambda <- par[1]; lambda.k <- par[2]
    c <- m * lambda / (m - 1)
    ck <- lambda.k - lambda / (m - 1)
    h <- 1 / (c * delta + 1)
    x1 <- x * sqrt(h)
    y1 <- y * sqrt(h)
    xx <- crossprod(x1)
    xy <- crossprod(x1, y1)
    yy <- crossprod(y1)
    yx <- t(xy)
    dd <- sum(log(1 / h))
    z1 <- z * sqrt(h)
    xz <- crossprod(x1, z1)
    yz <- crossprod(y1, z1)
    zz <- crossprod(z1)
    zx <- t(xz)
    zy <- t(yz)
    mm <- ck * solve(ck * zz + diag(q))
    xHx <- xx - xz %*% mm %*% zx
    xHy <- xy - xz %*% mm %*% zy
    xHz <- xz - xz %*% mm %*% zz
    yHz <- yz - yz %*% mm %*% zz
    zHz <- zz - zz %*% mm %*% zz
    yHy <- yy - yz %*% mm %*% zy
    yHx <- t(xHy)
    zHx <- t(xHz)
    zHy <- t(yHz)
    dd1 <- log(det(ck * zz + diag(q))) + dd  # Corrected on 12-11-2016
    yPy <- yHy - yHx %*% solve(xHx) %*% xHy
    dd2 <- log(det(xHx))
    loglike <- -0.5 * dd1 - 0.5 * dd2 - 0.5 * (n - r) * log(yPy)
    return(-loglike)
  }
  
  fixed <- function(par) {
    lambda <- par[1]
    lambda.k <- par[2]
    c <- m * lambda / (m - 1)
    ck <- lambda.k - lambda / (m - 1)
    h <- 1 / (c * delta + 1)
    x1 <- x * sqrt(h)
    y1 <- y * sqrt(h)
    xx <- crossprod(x1)
    xy <- crossprod(x1, y1)
    yy <- crossprod(y1)
    yx <- t(xy)
    dd <- sum(log(1 / h))
    z1 <- z * sqrt(h)
    xz <- crossprod(x1, z1)
    yz <- crossprod(y1, z1)
    zz <- crossprod(z1)
    zx <- t(xz)
    zy <- t(yz)
    mm <- ck * solve(ck * zz + diag(q))
    xHx <- xx - xz %*% mm %*% zx
    xHy <- xy - xz %*% mm %*% zy
    xHz <- xz - xz %*% mm %*% zz
    yHz <- yz - yz %*% mm %*% zz
    zHz <- zz - zz %*% mm %*% zz
    yHy <- yy - yz %*% mm %*% zy
    yHx <- t(xHy)
    zHx <- t(xHz)
    zHy <- t(yHz)
    yPy <- yHy - yHx %*% solve(xHx) %*% xHy
    s2 <- drop(yPy / (n - r))
    beta <- solve(xHx, xHy)
    gamma <- lambda.k * (zHy - zHx %*% beta)
    v <- (lambda.k * diag(q) - lambda.k * zHz * lambda.k) * s2
    stderr <- sqrt(diag(v))
    result <- list(beta, gamma, stderr, s2, v)
    return(result)
  }
  
  par0 <- c(1, 1)
  #parm<-optim(par=par0,fn=func,method="Brent",lower=1e-8,upper=1e8)
  parm <- optim(par = par0, fn = func, method = "L-BFGS-B", lower = 1e-8, upper = 1e8)
  par <- parm$par
  effect <- fixed(par = par)
  beta <- effect[[1]]; gamma <- effect[[2]]; stderr <- effect[[3]]
  s2 <- effect[[4]]; cov <- effect[[5]]
  covparm <- c(par * s2, s2, par)
  
  l1 <- -func(par)
  par0 <- c(1, 0)
  parm0 <- optim(par = par0, fn = func, method = "L-BFGS-B", lower = c(1e-8, 1e-9), upper = c(1e8, 1e-8))
  par0 <- parm0$par
  l0 <- -func(par0)
  lrt <- -2 * (l0 - l1)
  
  out <- list(covparm, beta, gamma, stderr, cov, lrt)
  return(out)
}


random1 <- function(x, z, y, m, kk = qq, cov = "qq", delta, n) {
  if (cov == "kk") {
    qq <- eigen(kk, symmetric = T)
    uu <- qq$vectors; delta <- qq$values
    x <- t(uu) %*% x; y <- t(uu) %*% y; z <- t(uu) %*% z
  }
  r <- ncol(x)
  q <- ncol(z) / 3
  indx <- (q + 1):(q + q)
  func <- function(par) {
    lambda <- par[1]
    lambda.k <- par[2]
    c <- m * lambda / (m - 3)
    c1 <- rep(-lambda / (m - 3), q)
    ck <- rep(lambda.k - lambda / (m - 3), q)
    Ck <- diag(c(c1, ck, c1))
    h <- 1 / (c * delta + 1)
    x1 <- x * sqrt(h)
    y1 <- y * sqrt(h)
    xx <- crossprod(x1)
    xy <- crossprod(x1, y1)
    yy <- crossprod(y1)
    yx <- t(xy)
    dd <- sum(log(1 / h))
    z1 <- z * sqrt(h)
    xz <- crossprod(x1, z1)
    yz <- crossprod(y1, z1)
    zz <- crossprod(z1)
    zx <- t(xz)
    zy <- t(yz)
    mm <- Ck %*% solve(zz %*% Ck + diag(3 * q))
    xHx <- xx - xz %*% mm %*% zx
    xHy <- xy - xz %*% mm %*% zy
    xHz <- xz - xz %*% mm %*% zz
    yHz <- yz - yz %*% mm %*% zz
    zHz <- zz - zz %*% mm %*% zz
    yHy <- yy - yz %*% mm %*% zy
    yHx <- t(xHy)
    zHx <- t(xHz)
    zHy <- t(yHz)
    dd1 <- log(det(zz %*% Ck + diag(3 * q))) + dd
    yPy <- yHy - yHx %*% solve(xHx) %*% xHy
    dd2 <- log(det(xHx))
    loglike <- -0.5 * dd1 - 0.5 * dd2 - 0.5 * (n - r) * log(yPy)
    return(-loglike)
  }
  
  fixed <- function(par) {
    lambda <- par[1]
    lambda.k <- par[2]
    c <- m * lambda / (m - 3)
    c1 <- rep(-lambda / (m - 3), q)
    ck <- rep(lambda.k - lambda / (m - 3), q)
    Ck <- diag(c(c1, ck, c1))
    h <- 1 / (c * delta + 1)
    x1 <- x * sqrt(h)
    y1 <- y * sqrt(h)
    xx <- crossprod(x1)
    xy <- crossprod(x1, y1)
    yy <- crossprod(y1)
    yx <- t(xy)
    dd <- sum(log(1 / h))
    z1 <- z * sqrt(h)
    xz <- crossprod(x1, z1)
    yz <- crossprod(y1, z1)
    zz <- crossprod(z1)
    zx <- t(xz)
    zy <- t(yz)
    mm <- Ck %*% solve(zz %*% Ck + diag(3 * q))
    xHx <- xx - xz %*% mm %*% zx
    xHy <- xy - xz %*% mm %*% zy
    xHz <- xz - xz %*% mm %*% zz
    yHz <- yz - yz %*% mm %*% zz
    zHz <- zz - zz %*% mm %*% zz
    yHy <- yy - yz %*% mm %*% zy
    yHx <- t(xHy)
    zHx <- t(xHz)
    zHy <- t(yHz)
    yPy <- yHy - yHx %*% solve(xHx) %*% xHy
    s2 <- drop(yPy / (n - r))
    beta <- solve(xHx, xHy)
    gamma <- lambda.k * (zHy[indx, ] - zHx[indx, ] %*% beta)
    v <- (lambda.k * diag(q) - lambda.k * zHz[indx, indx] * lambda.k) * s2
    stderr <- sqrt(diag(v))
    result <- list(beta, gamma, stderr, s2, v)
    return(result)
  }
  
  par0 <- c(1, 1)
  #parm<-optim(par=par0,fn=func,method="Brent",lower=1e-8,upper=1e8)
  parm <- optim(par = par0, fn = func, method = "L-BFGS-B", lower = 1e-8, upper = 1e8)
  par <- parm$par
  effect <- fixed(par = par)
  beta <- effect[[1]]
  gamma <- effect[[2]]
  stderr <- effect[[3]]
  s2 <- effect[[4]]
  cov <- effect[[5]]
  covparm <- c(par * s2, s2, par)
  
  l1 <- -func(par)
  par0 <- c(1, 0)
  parm0 <- optim(par = par0, fn = func, method = "L-BFGS-B", lower = c(1e-8, 1e-9), upper = c(1e8, 1e-8) )
  par0 <- parm0$par
  l0 <- -func(par0)
  lrt <- -2 * (l0 - l1)
  
  out <- list(covparm, beta, gamma, stderr, cov, lrt)
  return(out)
}
