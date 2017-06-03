\name{plotBinmap}
\alias{plotBinmap}

\title{Make a plot of Binmap}

\usage{
   plotBinmap(genotype="", cols=c("red", "blue", "grey40"), ...)
}

\arguments{
    \item{genotype}{A data frame representing the genotype of a population.}
    \item{cols}{Colours used to indicate different genotypes.}
    \item{...}{parameters passed to plot.}
}

\description{
  Each row of the genotype data frame represents a Bin. The first column is the 
  ID of each Bin. The second column is the chromosome ID of each Bin. The 3rd 
  column is the start coordinate of each Bin. The 4th column is the end coordinate 
  of each Bin. Each of the rest columns gives the genotype of each accession 
  at different bins. The names of the rest columns can be any characters 
  specified by the User.
}

\author{Wen Yao}

\examples{
 
 ril.geno <- read.csv(system.file("examples/ril.geno.csv", package="binQTL"), as.is=TRUE)
 dim(ril.geno)
 ril.geno[1:2, 1:9]
 ## Not run
 # par(mar = c(3, 4, 1, 1))
 # plotBinmap(ril.geno, xlab="Genomic position", ylab="RILs", cex.axis=0.6)
 ## Not run
 
}
