\name{coupling.fields}
\alias{coupling.fields}
\title{Calculate coupling fields for DCA}
\description{
Calculates coupling fields for DCA from the connected-correlation matrix.
}
\usage{coupling.fields(C,out.aa="-",naa=21,npos=ncol(C)/naa)
}
\arguments{
  \item{C}{connected-correlation matrix as the one returned by Cmat}
  \item{out.aa}{AminoAcid taken as referenced, and removed from initial calculations. See references for more information.}
  \item{naa}{Number of aminoacids present on the MSA, default 21.}
  \item{npos}{Length of sequences at the MSA (ie, number of columns at the MSA).}
}
\value{
List with the following components:
\item{$eij}{matrix with coupling energies.}
\item{$eij.balanced}{matrix with coupling energies, such that sum_aa e_ij =1}
\item{$Cdiag}{Some technichal stuff on the calculations.}
}
\references{
Morcos, Faruck, et al. "Direct-coupling analysis of residue coevolution captures native contacts across many protein families." Proceedings of the National Academy of Sciences 108.49 (2011): E1293-E1301.

Espada et al.
}
\seealso{Cmat}
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

  ## Calculate coupling fields
  cf<-coupling.fields(Cmatrix)
  
  cf$eij
  cf$eij.balanced
}
\keyword{DCA}\keyword{local}\keyword{fields}
