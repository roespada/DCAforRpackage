\name{freq.joint.correction}
\alias{freq.joint.correction}
\title{Correction for joint frequency in weighted MSA}
\description{
  Corrects joint frequencies of aminoacids in each position (column) of a multiple sequence alignment (MSA).
}
\usage{
  freq.joint.correction(f,Mef,naa=21,npos=ncol(f)/naa,lambda=Mef)
}
\arguments{
  \item{f}{Matrix with joint frequencies of residues on each position of the MSA, as the one returned by aa.freq.joint.}
  \item{Mef}{Number of effective sequences on the original MSA. Recommended values: sum the weights given to each sequence, or value returned by Henikoff.w$w. }
  \item{npos}{Length of sequences at the original MSA. Default ncol(f)/naa.}
  \item{lambda}{Pseudocounter used to correct frequencies. Default is equivalent to Mef. }
  \item{naa}{number of aminoacids used. Default 21.}
}
\value{
  Matrix with corrected joint frequencies between pairs of aminoacids-at-positions. Columns/rows goes over positions, and each position over aminoacids. Ie, they go pos1.aa1, pos1.aa1, ... , pos1.aa21, pos2.aa1, pos2.aa2, ...  See example below.
    The correction is made according to 
    \deqn{latex}{f.corrected=((\lambda/naa^2)+f)/(\lambda+Mef) if pos i != pos j}
    \deqn{latex}{f.corrected=f/(\lambda+Mef) if pos i = pos j and aa i != aa j}
    \deqn{latex}{f.corrected=(f+(\lambda/naa))/(\lambda+Mef) if pos i = pos j and aa i = aa j}
}
\references{
Morcos, Faruck, et al. "Direct-coupling analysis of residue coevolution captures native contacts across many protein families." Proceedings of the National Academy of Sciences 108.49 (2011): E1293-E1301.
}
\seealso{
  \code{\link{aa.freq.joint},\link{freq.marg.correction}}
}
\examples{
  ##Example of MSA
  msa.letters<-rbind(c("G","Y","V","G","S"),c("G","F","D","G","F"),c("G","Y","D","G","F"),c("G","Y","Q","G","G"))
  msa<-msa2num(msa.letters)
    
  ##Calculate joint frequencies, sequences weighted according Henikoff:
  fj<-aa.freq.joint(msa,weight=Henikoff.w(msa)$w)
  
  ## Corrected joint frequencies
  freq.joint.correction(fj,Mef=Henikoff.w(msa)$mef)
}
\keyword{frequency}
