\name{readAlignment}
\alias{readAlignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Read a multiple sequence alignment in fasta format.
}
\description{
Read a multiple sequence alignment in fasta format and returns a matrix in 1-letter code.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
readAlignment(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{File in fasta format to read.
%%     ~~Describe \code{file} here~~
}
}

\examples{

file<-system.file("PF00023_ncbi.txt", package = "DCAforR")
msa<-readAlignment(file)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ alignment }
\keyword{ msa }
