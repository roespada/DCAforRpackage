# DCAforRpackage
R package for mfDCA and mfDCAid calculations


# Installation

You must have R installed.
If you do not have devtools installed, you can get it from CRAN.
In R type:

     library(devtools)

     install_github("roespada/DCAforRpackage",subdir="DCAforR")

     library(DCAforR)

# Example:

\# Guide to calculate DIid for a multiple sequence alignment from PFAM

require(DCA)

\# Load alignment and prepare the 2-repeats alignment

\# Please notice that a txt file will be saved on current working directory

msa<-pfam("PF00023",alignment = "ncbi")

MSA.n.neighbours(msa = msa,n = 1,outfile = "PF00023_2reps.txt")

\# DI and DIid calculations - this can take several minutes (10 to 15 min each DI)

DIid<-CompleteDIrepeat(msa = "PF00023_2reps.txt",cutoff.gaps = 0.7)

\# Plot DI and DIid matrix

image(DIid$di,x=1:nrow(DIid$di),y=1:ncol(DIid$di),col=grey(seq(1,0,by=-0.001)),
      xlab="Position",ylab="Position",main="DIid");box()





 
 
