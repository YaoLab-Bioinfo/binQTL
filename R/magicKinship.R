
######################### calculated kinship
magicKinship <- function(geno="") {
  gen <- geno[, -c(1:4)]
  m <- nrow(gen)
  n <- ncol(gen)
  
  kinship <- matrix(0, n, n)
  for(i in 1:(m/8)){
    dummy <- as.matrix(gen[(8*i-7):(8*i), ])
    kin <- crossprod(dummy)
    kinship <- kinship+kin
  }
  
  kinship.stander <- kinship/mean(diag(kinship))
  return(kinship.stander)
}
