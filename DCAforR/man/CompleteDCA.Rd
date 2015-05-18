\name{CompleteDCA}
\alias{CompleteDI}
\alias{CompleteDIid}
\title{Complete DCA calculations over a MSA}
\description{Calculates Direct Information from the marginal frequencies and coupling fields.}
\usage{
CompleteDI(msa,cutoff.gaps=0.7)
CompleteDIid(msa,cutoff.gaps=0.7)
}
\arguments{
  \item{msa}{Name of the fasta file.}
  \item{cutoff.gaps}{Cut-off for maximum number of gaps allowed on a position of the MSA.}
}
\value{
A matrix with DI or DIid values for each pair of positions.
}
\references{
Espada, R., Parra, R. G., Mora, T., Walczak, A. M., & Ferreiro, D. (2014). Capturing coevolutionary signals in repeat proteins. arXiv preprint arXiv:1407.6903.

Morcos, F., Pagnani, A., Lunt, B., Bertolino, A., Marks, D. S., Sander, C., ... & Weigt, M. (2011). Direct-coupling analysis of residue coevolution captures native contacts across many protein families. Proceedings of the National Academy of Sciences, 108(49), E1293-E1301.
}

\examples{
require(DCA)

#Load alignment and prepare the 2-repeats alignment
# Please notice that a txt file will be saved on the current working directory
msa<-pfam("PF00023",alignment = "ncbi")
MSA.n.neighbours(msa = msa,n = 1,outfile = "PF00023_2reps.txt")


#DI and DIid calculations - this can take several minutes (10 to 15 min each DI)
DI<-CompleteDI(msa = "PF00023_2reps.txt",cutoff.gaps = 0.7)
DIid<-CompleteDIrepeat(msa = "PF00023_2reps.txt",cutoff.gaps = 0.7)


#Plot DI and DIid matrix
image(DI$DI,x=1:nrow(DI$DI),y=1:ncol(DI$DI),col=grey(seq(1,0,by=-0.001)),
      xlab="Position",ylab="Position",main="DI");box()
image(DIid$di,x=1:nrow(DIid$di),y=1:ncol(DIid$di),col=grey(seq(1,0,by=-0.001)),
      xlab="Position",ylab="Position",main="DIid");box()



}
\keyword{di}\keyword{direct}\keyword{information}\keyword{dca}\keyword{mutual}
