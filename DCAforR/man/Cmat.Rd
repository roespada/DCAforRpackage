\name{Cmat}
\alias{Cmat}
\title{Connected-correlation matrix}
\description{Calculates connected-correlation matrix from the matrices of marginal frequencies and joint frequencies.
}
\usage{Cmat(fmar,fjoint)
}
\arguments{
  \item{fmar}{Matrix of marginal frequencies, as the one returned by aa.freq.marg or freq.marg.correction}
   \item{fjoint}{Matrix of joint frequencies, as the one returned by aa.freq.joint or freq.joint.correction}
}
\value{
 Matrix with connected-correlation matrix between pairs of aminoacids-at-positions. Columns/rows goes over positions, and each position over aminoacids. Ie, they go pos1.aa1, pos1.aa1, ... , pos1.aa21, pos2.aa1, pos2.aa2, ...  See example below.
  We define connected-correlation matrix as
  \deqn{latex}{C= f_{ij}(A_i,A_j)-f_i(A_i)f_j(A_j)}
  where \eqn{latex}{f_{ij}(A_i,A_j)} is the joint frequency of aminoacid \eqn{latex}{A_i} in position i and position \eqn{latex}{A_j} in position j, and \eqn{latex}{f_{i}(A_i)} is the marginal frequency of aminoacid \eqn{latex}{A_i} in positoin i See references for further explanation.
 
}
\references{
Morcos, Faruck, et al. "Direct-coupling analysis of residue coevolution captures native contacts across many protein families." Proceedings of the National Academy of Sciences 108.49 (2011): E1293-E1301.
}
\seealso{
  \code{}
}
\examples{
  ##Example of MSA
  msa.letters<-rbind(c("G","Y","V","G","S"),c("G","F","D","G","F"),c("G","Y","D","G","F"),c("G","Y","Q","G","G"))
  msa<-msa2num(msa.letters)
    
  ##Calculate marginal frequencies, sequences weighted according Henikoff:
  fm<-aa.freq.marg(msa,weight=Henikoff.w(msa)$w)
  ## Corrected marginal frequencies
  fmc<-freq.marg.correction(fm,Mef=Henikoff.w(msa)$mef)
  
  ##Calculate joint frequencies, sequences weighted according Henikoff:
  fj<-aa.freq.joint(msa,weight=Henikoff.w(msa)$w)
  ## Corrected joint frequencies
  fjc<-freq.joint.correction(fj,Mef=Henikoff.w(msa)$mef)
  
  ## Calculate connected-correlation matrix
  Cmatrix<-Cmat(fmc,fjc)
  Cmatrix
}
\keyword{}
