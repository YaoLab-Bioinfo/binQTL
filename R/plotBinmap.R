
plotBinmap <- function(genotype="", cols=c("red", "blue", "grey40"), ...) {
  chr.len <- tapply(genotype$Stop, genotype$Chr, max)
  x.brks <- cumsum(chr.len) - chr.len/2
  chr.cum.len <- c(0, cumsum(chr.len)[-length(chr.len)])
  names(chr.cum.len) <- names(chr.len)
  genotype$Start <- genotype$Start + chr.cum.len[genotype$Chr]
  genotype$Stop <- genotype$Stop + chr.cum.len[genotype$Chr]
  
  plot(5, type="n", xlim=c(0, max(genotype$Stop)),
       ylim=c(0, ncol(genotype)-4), xaxt="n", ...)
  axis(side=1, at=x.brks, labels=unique(genotype$Chr))
  
  names(cols) <- c("P1", "P2", "H")
  
  for (i in 5:ncol(genotype)) {
    rect(xleft=genotype$Start, ybottom = i-5, xright=genotype$Stop,
         ytop=i-4, col=cols[genotype[, i]], border=NA)
  }
  
  abline(v=c(0, cumsum(chr.len)), col="grey80", lwd=0.5)
}
