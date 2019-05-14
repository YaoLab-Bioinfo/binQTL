
plotQTL <- function(qtl.res="", y.bottom=0, y.top=2, cols=c("grey30", "grey70"), ylab="-logp", xlab="", ...) {
  qtl.res$Start <- as.numeric(qtl.res$Start)
  qtl.res$Stop <- as.numeric(qtl.res$Stop)
  chr.len <- tapply(qtl.res$Stop, qtl.res$Chr, max)
  x.brks <- cumsum(chr.len) - chr.len/2
  chr.cum.len <- c(0, cumsum(chr.len)[-length(chr.len)])
  names(chr.cum.len) <- names(chr.len)
  qtl.res$Start <- qtl.res$Start + chr.cum.len[qtl.res$Chr]
  qtl.res$Stop <- qtl.res$Stop + chr.cum.len[qtl.res$Chr]
  qtl.res$col <- cols[1]
  if (length(chr.len)>=2) {
    qtl.res$col[qtl.res$Chr %in% names(chr.len)[seq(2, length(chr.len), by=2)]] <- cols[2]
  }

  plot(5, type="n", xlim=c(0, max(qtl.res$Stop)),
       ylim=c(y.bottom, max(qtl.res[,5])+y.top), ylab=ylab, xlab=xlab, xaxt="n", ...)
  axis(side=1, at=x.brks, labels=unique(qtl.res$Chr))
  rect(xleft=qtl.res$Start, ybottom = 0, xright=qtl.res$Stop,
       ytop=qtl.res[,5], col=qtl.res$col, border=NA)
}
