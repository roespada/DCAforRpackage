\name{DCAforR-package}
\alias{DCAforR-package}
\alias{DCAforR}
\docType{package}
\title{
Calculation of Direct Coupling Analysis (mfDCA) for repeat proteins.
}
\description{
This package is an implementation of mfDCA and mfDCAid in R. You can reproduce the results in the paper cited below with it.

}
\details{
\tabular{ll}{
Package: \tab DCAforR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-13\cr
License: \tab MIT + file LICENSE\cr
}

This package includes functions to load a MSA in fasta format, translate it to numerical code, generate alignment of pairs of repeats (resticted to PFAM ncbi alignments), calculate marginal and joint frequencies of residues, weight sequences according to Henikoff weights, and calculate DI and DIid. The easiest way to use it is to run directly \code{CompleteDI} and \code{CompleteDIid}
}
\author{Rocio Espada

Maintainer: Rocio Espada <roespada@qb.fcen.uba.ar>
}
\references{
Espada, R., Parra, R. G., Mora, T., Walczak, A. M., & Ferreiro, D. (2015). Capturing coevolutionary signals in repeat proteins. arXiv preprint arXiv:1407.6903.
}
\keyword{ direct }
\keyword{ information }
\keyword{ protein }
\keyword{ repetitive protein }

\examples{

#Load alignment and prepare the 2-repeats alignment
# Please notice that a txt file will be saved on the current working directory
file<-system.file("PF00023_ncbi.txt", package = "DCAforR")
msa<-readAlignment(file)
MSA.n.neighbours(msa = msa,n = 1,outfile = "PF00023_2reps.txt")


#DI and DIid calculations - this can take several minutes (10 to 15 min each DI)
DI<-CompleteDI(msa = "PF00023_2reps.txt",cutoff.gaps = 0.7)
DIid<-CompleteDIid(msa = "PF00023_2reps.txt",cutoff.gaps = 0.7)


#Plot DI and DIid matrix
image(DI$DI,x=1:nrow(DI$DI),y=1:ncol(DI$DI),col=grey(seq(1,0,by=-0.001)),
      xlab="Position",ylab="Position",main="DI");box()
image(DIid$di,x=1:nrow(DIid$di),y=1:ncol(DIid$di),col=grey(seq(1,0,by=-0.001)),
      xlab="Position",ylab="Position",main="DIid");box()



}
