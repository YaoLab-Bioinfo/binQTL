\name{plotQTL}
\alias{plotQTL}

\title{Make a plot of QTL mapping result}

\usage{
   plotQTL(qtl.res="", y.bottom=0, y.top=2, cols=c("grey30", "grey70"),
           ylab="-logp", xlab="", ...)
}

\arguments{
    \item{qtl.res}{A data frame of QTL scanning result.}
    \item{y.bottom}{A numeric number indicating the min y coordinate range.}
    \item{y.top}{A numeric number used to alter the max y coordinate range.}
    \item{cols}{Colours used to indicate different chromosome.}
    \item{ylab}{Y axis label.}
    \item{xlab}{X axis label.}
    \item{...}{parameters passed to plot.}
}

\description{
  The qtl.res data frame should contain 5 columns. Each row of this data frame
  represents a QTL. The 1st column is the ID of a Bin representing a QTL. 
  The 2nd column is the chromosome ID of each QTL. The 3rd column is the start 
  coordinate of a QTL. The 4th column is the end coordinate of a QTL. The 
  5th column is the negative Log P value of QTL scanning.
}

\author{Wen Yao}

\examples{
 
  ril.phe <- read.csv(system.file("examples/ril.phe.csv", package="binQTL"), as.is=TRUE)
  dim(ril.phe)
  ril.geno <- read.csv(system.file("examples/ril.geno.csv", package="binQTL"), as.is=TRUE)
  dim(ril.geno)
  qtl.res <- aovQTL(phenotype = ril.phe, genotype = ril.geno)
  plotQTL(qtl.res[, c(1:4, 6)], main="ANOVA", cex.axis=0.6, ylab=expression(-log[10](p)))
 
}
