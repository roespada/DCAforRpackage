###############################################################################################
# General functions
###############################################################################################

#Load Alignment
#file: name of the fasta with the MSA
#returns: matrix with alignment
readAlignment<-function(file){
  lineas<-readLines(file,warn = FALSE)
  lineas.nombres<-grep(pattern = ">",x=lineas)
  
  nc <- sum(sapply(lineas[(lineas.nombres[1]+1) : (lineas.nombres[2]-1) ],nchar))
  msa<-matrix(nrow=length(lineas.nombres),ncol=nc)
  
  for(i in seq_along(lineas.nombres)){
    lin.seq.init<-lineas.nombres[i]+1
    if(i == length(lineas.nombres)){
      lin.seq.end<-length(lineas)
    }else{
      lin.seq.end<-lineas.nombres[i+1]-1
    }
    
    complete.seq.protein <- paste(lineas[lin.seq.init:lin.seq.end],collapse = "")
    complete.seq.protein <- unlist(strsplit(complete.seq.protein,split = ""))
    
    msa[i,]<-complete.seq.protein
  }
    rownames(msa)<-sub(pattern=">",replacement = "",lineas[lineas.nombres])
  
  return(msa)
}

## Alignment msa to integer numbers matrix
#msa: msa matrix as given by read.fasta()$ali
msa2num<-function(msa){
  if(is.matrix(msa)){ids<-rownames(msa)}
  aminoac<-c("-","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  msa<-toupper(msa)
  if(any(! msa%in%aminoac)){
    pos<-which(! msa%in%aminoac)
    warning(length(pos)," Positions had unidentified residues (",unique(msa[pos]),")\n \t All were set as gap")
    msa[pos]<-"-"
    
  }
  
  msa[which(msa=="-" | msa=="X"| msa==".")]<-1
  msa[which(msa=="A" | msa=="a")]<-2
  msa[which(msa=="C" | msa=="c")]<-3
  msa[which(msa=="D" | msa=="d")]<-4
  msa[which(msa=="E" | msa=="e")]<-5
  msa[which(msa=="F" | msa=="f")]<-6
  msa[which(msa=="G" | msa=="g")]<-7
  msa[which(msa=="H" | msa=="h")]<-8
  msa[which(msa=="I" | msa=="i")]<-9
  msa[which(msa=="K" | msa=="k")]<-10
  msa[which(msa=="L" | msa=="l")]<-11
  msa[which(msa=="M" | msa=="m")]<-12
  msa[which(msa=="N" | msa=="n")]<-13
  msa[which(msa=="P" | msa=="p")]<-14
  msa[which(msa=="Q" | msa=="q")]<-15
  msa[which(msa=="R" | msa=="r")]<-16
  msa[which(msa=="S" | msa=="s")]<-17
  msa[which(msa=="T" | msa=="t")]<-18
  msa[which(msa=="V" | msa=="v")]<-19
  msa[which(msa=="W" | msa=="w")]<-20
  msa[which(msa=="Y" | msa=="y")]<-21
  
  
  if(is.matrix(msa)){
    msa<-apply(msa,2,as.integer)
    rownames(msa)<-ids
  }else{msa<-as.integer(msa)}
  return(msa)
}


## Calculate Henikoff's weights
# # Example of the paper to check with...
# msa.o<-rbind(c("G","Y","V","G","S"),c("G","F","D","G","F"),c("G","Y","D","G","F"),c("G","Y","Q","G","G"))
# msa<-msa2num(msa.o)
#msa: integer matrix as out of msa2num to calculate the weights to.
Henikoff.w<-function(msa,fmarg,normalize=FALSE){
  #Renombro fmarg
  tabla<-fmarg
  
  #number of amino acids on each position (r)
  r<-sapply(apply(msa,2,unique),length)
  #number of sequences with each amino acid on each position on each sequnce (s)
  s<-sapply(1:ncol(msa),function(i,msa,vv){
    vv[i,msa[,i]]
  },msa=msa,vv=tabla)
  
  #calculate matrix of 1/rs
  w<-matrix(rep(1,length(msa)),ncol=ncol(msa),nrow=nrow(msa))
  w<-(t(w)/(t(s)*r))
  
  #sum to get the sequence weight
  w<-apply(w,2,sum)
  #normalize and calculate Mef
  if(normalize==T){
    w<-w/sum(w)
    mef<-exp(sum(-w*log(w)))
  }else{
    mef<-sum(w)
  }
  names(w)<-NULL
  if(!is.null(rownames(msa))){names(w)<-rownames(msa)}
  return(list(w=w,mef=mef))
}


## Remove positions with many gaps
msa.removed.gaps<-function(msa,cutoff.gaps=0.7,aa2rm="-"){
  gaps<-apply(msa,2,function(x){length(which(x==aa2rm))})/nrow(msa)
  msa<-msa[,which(gaps<cutoff.gaps)]
  return(msa)
}


###############################################################################################
# Calculation of frequencies
###############################################################################################
##  Marginal frequencies of amino acid in each position     NEEDS DESCR PACKAGE Y TIENE ALGO MEDIO RARO CON LOS AA...
#M: multiple sequence alignment - as integer matrix
#weights: vector of weights assigned to each sequence
aa.freq.marg<-function(M,weight=rep(1,nrow(M)),npos=dim(M)[2],naa=21,colnames.aa=T){
  # Verification:
  if(nrow(M)!=length(weight)) stop("Different number of sequences and weights!")
  
  #     if(is.null(npos)){npos<-dim(M)[2]
  #                       warning("npos set as dim(M)[2]")}
  #     if(is.null(naa)){naa<-length(unique(c(M)))
  #                      warning("naa set as length(unique(c(M)))")}
  
  
  frecuencias.marginales<-matrix(ncol=naa,nrow=npos)
  
  for(i in 1:npos){
    #tomo la posicion i esima del alineamiento
    aux<-freq(x=M[,i],w=weight,plot="FALSE")          		#calculo la frecuencia de cada aminoacido pesada por wieghts
    amP<-attributes(aux)$xdata.c    					    #me quedo el numero total
    frecuencias.marginales[i,as.numeric(names(amP))]<-amP				#lo meto en la matriz
  }
  frecuencias.marginales[which(is.na(frecuencias.marginales))]<-0					#si quedo algun NA lo mando a cero
  if(colnames.aa==T){
    colnames(frecuencias.marginales)<-c("-","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  }
  
  return(frecuencias.marginales)
}

##  Joint Frequencies of amino acids on two positions 
#M: multiple sequence alignment - as integer matrix
#weights: vector of weights assigned to each sequence
aa.freq.joint<-function(M,weight=rep(1,nrow(M)),naa=21,npos=dim(M)[2],colnames.aa=T){
  # Verification:
  if(nrow(M)!=length(weight)) stop("Different number of sequences and weights!")
  
  #     #Definiciones basicas del alineamiento
  #     if(is.null(npos)){npos<-dim(M)[2]
  #                       warning("npos set as dim(M)[2]")}
  #     if(is.null(naa)){naa<-length(unique(c(M)))
  #                      warning("naa set as length(unique(c(M)))")}
  #     
  frecuencias.conjuntas2<-matrix(nrow=npos*naa,ncol=npos*naa)
  for(i in 1:npos){
    ipos<-(naa*(i-1)+1):(naa*i)
    for(j in i:npos){
      aux<-xtabs(weight~M[,i]+M[,j])
      frecuencias.conjuntas2[ipos,(naa*(j-1)+1):(naa*j)][as.numeric(rownames(aux)),as.numeric(colnames(aux))]<-aux
    }
  }
  ft <- t(frecuencias.conjuntas2)
  ind <- lower.tri(frecuencias.conjuntas2,diag=FALSE)
  frecuencias.conjuntas2[ind] <- ft[ind]
  
  frecuencias.conjuntas2[which(is.na(frecuencias.conjuntas2))]<-0
  
  
  if(colnames.aa==T){
    aas<-c("-","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    n<-paste(rep(as.character(1:npos),each=naa),aas,sep=".")
    colnames(frecuencias.conjuntas2)<-rownames(frecuencias.conjuntas2)<-n
  }
  
  
  return(frecuencias.conjuntas2)
}

## Correction on marginal frequencies by pseudocounter and sequences weights
#f: marginal frequencies matrix. Columns:aminoacids Rows: Position on msa
#lambda: Correction factor
#Mef: sum of sequences' weights
freq.marg.correction<-function(f,Mef,lambda=Mef,naa=21){
  if(is.null(naa)){naa<-ncol(f)}
  out<-((lambda/naa)+f)/(lambda+Mef)
  return(out)
}

## Correction on joint frequencies by pseudocounter and sequences weights
#f: joint frequencies (matrix) 
#naa: number of possible amino acids (21 by default)
#npos: length of the sequences on the msa
#lambda: Correction factor
#Mef: sum of sequences' weights
freq.joint.correction<-function(f,Mef,naa=21,npos=ncol(f)/naa,lambda=Mef){
  
  out<-((lambda/naa^2)+f)/(lambda+Mef)
  for(i in 1:npos){   #esta es una correccion extra que le tuve que agregar a la diagonal de frecuencias.conjuntas.corregida para poder invertir
    aux.pos<-(naa*(i-1)+1):(naa*i)
    faux<-f[aux.pos,aux.pos]/(lambda+Mef)
    diag(faux)<-(diag(f[aux.pos,aux.pos])+(lambda/naa))/(lambda+Mef)
    out[aux.pos,aux.pos]<-faux
  }
  return(out)
}

## C matrix: fij-fi fj 
#fmar: corrected marginal frequencies (matrix)
#fconj: corrected joint frequencies (matrix)
Cmat<-function(fmar,fjoint){
  ##  Matriz C  ##
  v<-c(t(fmar))
  C<-fjoint-(v%*%t(v))
  
  return(C)
}



###############################################################################################
# Calculation of DI and MI
###############################################################################################
## Entropy per position
#freq.marg: marginal frequency
#weights: weights used for each sequence
#norm: if freq.marg should be normalized by the sym of weights.
Spos<-function(freq.marg,npos=nrow(freq.marg)){
  S<-apply(freq.marg,1,function(x){
    aux<-x*log2(x)
    aux[is.nan(aux)]<-0
    return(-sum(aux))
  })
  return(S)
}

##Arma una matriz random
#M: matriz de secuencias en las filas
# ex function.random.matrix
msa.indep<-function(M){
  Mrandom<-apply(M,2,function(x){x[sample.int(n=length(x))]})
  rownames(Mrandom)<-paste("SeqR",c(1:nrow(M)),sep=".")
  return(Mrandom)
}

## Invertir la matriz C y balancearla. Devuelve: inversa de C, diagonal de la inversa de C, Cinv balanceada
#C matriz a invertir
#out.aa: aminoacido que no uso en la inversion (aunque luego lo agrego como ceros en la matrix final)
coupling.fields<-function(C,out.aa="-",naa=21,npos=ncol(C)/naa){
  
  if(is.character(out.aa)){out.aa<-msa2num(out.aa)}
  
  outpos<-naa*c(0:(npos-1))+out.aa
  pos.sin.aa<-setdiff(1:nrow(C),outpos)
  Caux<-C[pos.sin.aa,pos.sin.aa]
  Caux2<-solve(Caux)
  
  Cinv<-matrix(rep(0,ncol(C)*nrow(C)),ncol=ncol(C),nrow=nrow(C))
  Cinv[pos.sin.aa,pos.sin.aa]<-Caux2
  Cinv<-Cinv*(-1)
  
  auxM<-matrix(rep(0,ncol(Cinv)*nrow(Cinv)),ncol=ncol(Cinv),nrow=nrow(Cinv))
  for(i in 1:npos){
    ipos<-(naa*(i-1)+1):(naa*(i))
    auxM[ipos,ipos]<-Cinv[ipos,ipos]
  }
  
  
  ## Balanceo la matriz de couplings Cinv (armo K)  ##
  K<-matrix(ncol=ncol(Cinv),nrow=nrow(Cinv))
  Ch<-(Cinv-auxM)
  for(i in 1:npos){
    indexi<-(naa*(i-1)+1):(naa*(i))
    for(j in 1:npos){
      indexj<-(naa*(j-1)+1):(naa*(j))
      Caux<-Ch[indexi,indexj]
      if(any(dim(Caux)==c(naa,naa))==FALSE){cat(i,"dim problems \n")}
      ja<-matrix(rep(apply(Caux,2,sum),nrow(Caux)),ncol=(naa),byrow=T)
      jb<-matrix(rep(apply(Caux,1,sum),ncol(Caux)),nrow=(naa),byrow=F)
      j<-sum(Caux)
      K[indexi,indexj]<-Caux-1/naa*ja-1/naa*jb+1/naa^2*j
    } 
  }
  
  if(length(colnames(C))>0){
    colnames(Cinv)<-rownames(Cinv)<-colnames(C)
    colnames(auxM)<-rownames(auxM)<-colnames(C)
    colnames(K)<-rownames(K)<-colnames(C)
  }
  return(list(eij=Cinv,Cdiag=auxM,eij.balanced=K))
}

## Direct Information calculations
#fmarg:frecuencias marginales corregidas
#K:campos de interaccion balanceados ( $eij.balanced de coupling fields)
DIfunc<-function(fmarg,K,naa=21,npos=nrow(fmarg)){
  
  fm.vec<-c(t(fmarg))
  
  #PROBABILIDADES "DIRECTAS"
  p.dir<-matrix(ncol=ncol(K),nrow=nrow(K))
  eps<-10
  ee<-0.5
  for(i in 1:(npos-1)){
    
    #ii<-colnames(K)[which(partir(colnames(K),keep=1)==i)]
    ii<-((naa)*(i-1)+1):((naa)*i)
    hi<-rep(0,naa)
    names(hi)<-ii
    
    for(j in (i+1):npos){
      #cat("           ",i,j,"\n",sep=" - ")
      #jj<-colnames(K)[which(partir(rownames(K),keep=1)==j)]
      jj<-((naa)*(j-1)+1):((naa)*j)
      hj<-rep(0,naa)
      names(hj)<-jj
      
      eij<-K[ii,jj]
      fia<-fm.vec[ii]
      fjb<-fm.vec[jj]
      eps<-10
      t<-0
      while(eps>10^(-6)){  
        t<-t+1
        
        #corrijo los valores de h_i(a) y h_j(b) segun cuanto difieren las frec calculadas con las medidas
        hi<-hi+ee*(fm.vec[ii]-fia)
        hj<-hj+ee*(fm.vec[jj]-fjb)
        
        #Calculo las sumas parciales        
        him<-t(matrix(rep(hi, naa),ncol=naa,byrow=T))
        #                 rownames(him)<-names(hi)
        hjm<-(matrix(rep(hj, naa),ncol=naa,byrow=T))
        #                 colnames(hjm)<-names(hj)
        
        x<-exp(eij+hjm+him)
        p.dir.aux<-x/sum(c(x))
        
        fia<-apply(p.dir.aux,1,sum)
        fjb<-apply(p.dir.aux,2,sum)
        
        #estimo el error que estoy cometiendo entre las frecuencias medidas fi fj y las calculadas con estos campos
        eps<-max(abs(c(fm.vec[ii]-fia,fm.vec[jj]-fjb))) 
      }
      
      p.dir[ii,jj]<-p.dir.aux
    }
  }
  p.dir[which(is.na(p.dir))]<-0
  p.dir<-p.dir+t(p.dir)
  
  #CALCULO DIRECT INFORMATION
  
  DI<-matrix(ncol=npos,nrow=npos)
  for(i in 1:npos){
    ii<-((naa)*(i-1)+1):((naa)*i)
    fi.inv<-1/fm.vec[ii]
    for(j in 1:npos){
      jj<-((naa)*(j-1)+1):((naa)*j)
      fj.inv<-1/fm.vec[jj]
      aux<-p.dir[ii,jj]*log(t(t(fi.inv*p.dir[ii,jj])*fj.inv))
      DI[i,j]<-sum(aux)
    }
  }
  DI[cbind(1:nrow(DI),1:nrow(DI))]<-0
  
  return(list(DI=DI,Pdir=p.dir))
}


###############################################################################################
# Functions for repeat proteins
###############################################################################################

# Armar alineamientos de vecinos n-esimos

#Funcion para armar alineamiento a n vecinos con datos bajados de PFAM/NCBI (n=1 primeros vecinos, n=2 segundos veinos, etc)
#maxdist=ncol(msa)*(n*1.5-1) Maxima cantidad de aminoacidos en el medio permitidas como inserciones.
#        el default es la cantidad de repeats en el medio + medio repeat por cada union de repeats: [ n-1 + n ] * ncol(msa)
MSA.n.neighbours<-function(msa,n=2,outfile=paste("MSA_neighb_",n,".fasta",sep=""),maxdist=ncol(msa)*(n*1.5-1),mindist=ncol(msa)*(n-1)){
  
  #ids de los repeats en matriz nombrePROT/inicioREP/finalREP
  nom<-rownames(msa)
  nom<-matrix(unlist(sapply(sapply(nom,strsplit,"/"),strsplit,"-")),ncol=3,byrow=T)
  
  ##### Busco los repeats que esten separados por n repeats  #######
  #calculo la frecuencia de ocurrencia por nombre de proteina (ver cuantos repeats tiene)
  f<-table(nom[,1])
  #busco los nombres de las proteinas que tienen al menos n+1 repeats asociados
  n1<-names(which(f>n))
  #en pares guardo los nombres de repeats a distancia n
  for(i in seq_along(n1)){
    enp<-which(nom[,1]==n1[i])  #indices de donde aparece la proteina n1[i]
    
    nom.sub<-nom[enp,]
    nom.sub<-nom.sub[order(as.numeric(nom.sub[,2])),]
    
    for(j in 1:(nrow(nom.sub)-n)){
      #si la distancia en secuencia es mayor a maxdist o los repeats se superponen, entonces lo salteo
      seqdist<-as.numeric(nom.sub[j+n,2])-as.numeric(nom.sub[j,3])
      if(seqdist>maxdist | seqdist< mindist ){next}
      #sino uno los dos repeats y lo agrego al alineamiento en el outfile
      r1<-paste(nom.sub[j,1],"/",nom.sub[j,2],"-",nom.sub[j,3],sep="")
      r2<-paste(nom.sub[j+n,1],"/",nom.sub[j+n,2],"-",nom.sub[j+n,3],sep="")
      
      cat(paste(">",r1,"-",r2,"\n",sep=""),file=outfile,append=T)
      cat(paste(c(msa[r1,],msa[r2,]),collapse = ""),file=outfile,append=T)
      cat("\n",file=outfile,append=T)
    }    
  }
  return(0)
}


# Peso para repeat proteins
weight.idREP<-function(msa.real){
  npos1<-ncol(msa.real)/2
  msa.mezcla<-cbind(msa.real[,1:npos1],msa.real[sample(1:nrow(msa.real)),(npos1+1):(2*npos1)])
  
  dist.primeros<-c()
  for(i in 1:nrow(msa.real)){
    dist.primeros<-c(dist.primeros,length(which(msa.real[i,1:npos1]==msa.real[i,(npos1+1):(2*npos1)]))/npos1)
  }
  names(dist.primeros)<-rownames(msa.real)
  if(any(is.na(dist.primeros))){dist.primeros[which(is.na(dist.primeros))]<-0}
  
  ## PID vecinos random
  dist.random<-c()
  for(i in 1:nrow(msa.mezcla)){
    dist.random<-c(dist.random,length(which(msa.mezcla[i,1:npos1]==msa.mezcla[i,(npos1+1):(2*npos1)]))/npos1)
  }
  names(dist.random)<-rownames(msa.mezcla)
  if(any(is.na(dist.random))){dist.random[which(is.na(dist.random))]<-0}
  
  
  bks<-seq(0,1,length.out=npos1)
  hrr<-hist(dist.random,breaks=bks,plot=F)
  h1<-hist(dist.primeros,breaks=bks,plot=F)
  
  # Armo pesos interpolando los histogramas
  w.seq.pid<-rep(NA,nrow(msa.real))
  cts<-h1$breaks
  pspR<-hrr$counts
  pspF<-h1$counts
  pspR<-pspR/sum(pspR)
  pspF<-pspF/sum(pspF)
  for( i in 1:(length(cts)-1) ){
    w.seq.pid[which(dist.primeros>=cts[i] & dist.primeros<(cts[i+1]+h1$mids[1]))]<- (pspR[i]/pspF[i])
  }
  names(w.seq.pid)<-rownames(msa.real)
  
  return(w.seq.pid)
  
}




###############################################################################################
# Functions for complete DI calculations
###############################################################################################

CompleteDI<-function(msa,cutoff.gaps=0.7){
  msa <- readAlignment(msa)
  msa <- msa.removed.gaps(msa,cutoff.gaps=cutoff.gaps)
  msa <- msa2num(msa)
  w <- Henikoff.w(msa,fmarg=aa.freq.marg(M=msa))
  
  fm  <- aa.freq.marg(M=msa,weight=w$w,colnames.aa=T)
  fmc <- freq.marg.correction(f=fm,Mef=w$mef)
  
  fj  <- aa.freq.joint(M=msa,weight=w$w,colnames.aa=T)
  fjc <- freq.joint.correction(f=fj,Mef=w$mef)
  
  c.matrix <- Cmat(fmar=fmc,fjoint=fjc)
  cf <- coupling.fields(C=c.matrix)
  di <- DIfunc(fmarg=fmc,K=cf$eij)
  
  return(di)
}

CompleteDIid<-function(msa,cutoff.gaps=0.7){
  msa <- readAlignment(msa)
  msa <- msa.removed.gaps(msa,cutoff.gaps=cutoff.gaps)
  msa <- msa2num(msa)
  w   <- Henikoff.w(msa,fmarg=aa.freq.marg(M=msa))
  wid <- weight.idREP(msa)
  
  w$w <- w$w*wid
  
  fm  <- aa.freq.marg(M=msa,weight=w$w,colnames.aa=T)
  fmc <- freq.marg.correction(f=fm,Mef=sum(w$w))
  
  fj  <- aa.freq.joint(M=msa,weight=w$w,colnames.aa=T)
  fjc <- freq.joint.correction(f=fj,Mef=sum(w$w))
  
  c.matrix <- Cmat(fmar=fmc,fjoint=fjc)
  cf <- coupling.fields(C=c.matrix)
  di <- DIfunc(fmarg=fmc,K=cf$eij)
  
  return(di=di$DI)
}

