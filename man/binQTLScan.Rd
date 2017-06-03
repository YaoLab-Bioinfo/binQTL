\name{binQTLScan}
\alias{binQTLScan}

\title{QTL mapping with binmap data}

\usage{
   binQTLScan(phenotype="", genotype="", population="RIL")
}

\arguments{
    \item{phenotype}{A data frame representing the phenotype of a population.}
    \item{genotype}{A data frame representing the genotype of a population.}
    \item{population}{A character string indicating the type of the population.}
}

\description{
  The phenotype data frame should contain two columns. The first column is the ID 
  of each accession of a population. The second column contains the phenotypic 
  value of specific trait for each accession. 
  Each row of the genotype data frame represents a Bin. The first column is the 
  ID of each Bin. The second column is the chromosome ID of each Bin. The 3rd 
  column is the start coordinate of each Bin. The 4th column is the end coordinate 
  of each Bin. Each of the rest columns gives the genotype of each accession 
  at different bins. The names of the rest columns can be any characters 
  specified by the User.
  The population type can be RIL or F2.
}

\author{Wen Yao, Shizhong Xu}

\examples{

  ril.phe <- read.csv(system.file("examples/ril.phe.csv", package="binQTL"), as.is=TRUE)
  dim(ril.phe)
  head(ril.phe)
  ril.geno <- read.csv(system.file("examples/ril.geno.csv", package="binQTL"), as.is=TRUE)
  dim(ril.geno)
  ril.geno[1:2, 1:9]
  qtl.res <- binQTLScan(phenotype = ril.phe, genotype = ril.geno, population="RIL")
  head(qtl.res)
  
}
