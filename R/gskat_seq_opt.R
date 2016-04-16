# GSKAT optimal test

##null fit function: estimate aplha.hat and delta.hat
fit_FSKAT_IC<-function(y,XC,FID){
  #y:phenotype
  #XC:covaraite matrix, including PC and intercept !
  #FID: family ID, not used if indentity working correlation
  
  #alpha<-t(t(c(-1,rep(1,(ncol(XC)-1))))) #initial values
  alpha<-summary(glm(y~XC[,-1],family=binomial("logit")))$coeff[,1]
 
  diff1<-10
  i<-0
  while (diff1>10e-6)
    {
    
    i<-i+1
    #Fisher scoring
    mu<- plogis( XC %*% alpha)
    mu=as.vector(mu)
    #D=diag((mu)*(1-mu))%*%X
    V=diag((mu)*(1-mu)) #V=sqrt(A)%*%(R)%*%sqrt(A) 
    U=t(XC)%*%t(t((y-mu)))    #U=t(D)%*%solve(V)%*%t(t((y-mu)))
    I=t(XC)%*%V%*%XC   	#I<-t(D)%*%solve(V)%*%D
    
    alpha.new<-alpha+solve(I)%*%U
    diff1 <- max(abs(alpha.new-alpha))
    alpha<-alpha.new

   }	
  return(list(alpha,i))
}

gskat_seq_opt<-function(y,XC,Z,ID,NB,impute.method="fixed",SNP.weights=NULL,w_a=1,w_b=25,resampling=TRUE,pw="Rade",Uc=TRUE,sW=FALSE,np=10000) 
{
  maf<-apply(Z,2,function(x){sum(x,na.rm=T)/((length(x)-sum(is.na(x)))*2)} )#cacluate before center
  if(impute.method=="fixed") Z<-apply(Z,2,function(x){x[is.na(x)]=mean(x,na.rm=TRUE);return(x)}) #mean values imputation
  
  Z<-scale(Z,scale=F) #Z need to be centered !!!
  if (sum(XC[,1]!=1)!=0) {XC=cbind(1,data.matrix(XC))} #now with intercept, assume all nonvariants removed
  Z<-data.matrix(Z)  #Z is the combined genotype, must be matrix format
  y<-as.vector(y);FID<-as.factor(ID$FID)
  
  ##sort according to FID	
  FID.old<-FID
  FID<-FID[order(FID.old)]; y<-y[order(FID.old)];XC<-XC[order(FID.old),];Z<-data.matrix(Z[order(FID.old),]) #!!
  ID<-ID[order(FID.old),]
  p=ncol(Z); q=ncol(XC)-1
  #
  #weight W
  #w_a=1; w_b=25
  if (is.null(SNP.weights))
    {
    #w_a=1; w_b=25
    if (length(maf)==1){W=1} else {
      #W<-diag((dbeta(maf,w_a, w_b))^2) } #beta density
      W<-diag(dbeta(maf,w_a, w_b)) } #beta density
    if (sW==TRUE) {W<-W/max(W)} 
    #W<-diag(length(maf))
    } else 
    {
    W<-diag(SNP.weights)
    }
  #
  null<-fit_FSKAT_IC(y,XC,FID)
  alpha<-null[[1]]
  #
  mu<-as.vector(plogis(XC%*%alpha))
  U=t(Z)%*%t(t((y-mu)))
  V=diag((mu)*(1-mu))
  Covy<-tcrossprod(y-mu)
  #Covy[((row(Covy)-1)%/%4+1)!=((col(Covy)-1)%/%4+1)]=0  #for same family strut of 4 members
  Covy=Covy*blockMatrixDiagonal(lapply(split(FID,FID),function(alist){matrix(1,length(alist),length(alist))}))
  XZ<-cbind(XC,Z)
  B<-t(XZ)%*%Covy%*%XZ  #Bi<-t(Di)%*%iVi%*%Covyi%*%iVi%*%Di
  A<-t(XZ)%*%V%*%XZ  #AAi<-t(Di)%*%iVi%*%Di 
  
  Azx<-A[(q+1+1):ncol(A),1:(q+1)]
  Axx<-A[1:(q+1),1:(q+1)]
  C<-cbind(-Azx%*%solve(Axx),diag(p))
  
  rho=seq(0,1,0.1) ;TS=rho ;pval=rho ;pvalb=rho;
  onev=matrix(1,nrow(W),1)
  #TS<-t(U)%*%((1-rho)*W+rho*tcrossprod(sqrt(W)%*%onev))%*%U    #new ts under rho
  for(i in 1:11)
  {TS[i]=t(U)%*%((1-rho[i])*W+rho[i]*tcrossprod(sqrt(W)%*%onev))%*%U}
  Bsqrt=denman.beavers(B)$sqrt
  BC<-Bsqrt%*%t(C)%*%((1-rho[i])*W+rho[i]*tcrossprod(sqrt(W)%*%onev))%*%C%*%Bsqrt  #to get new eigenvalues
  Lamda<-eigen(BC,only.values=T)$values
  results<-davies(TS,Lamda,rep(1, length(Lamda)))
  results
  if (resampling==TRUE | results$Qq>1 | results$ifault==1)
  {
    a=diag((mu)*(1-mu))
    R<-diag(dim(a)[1])
    iV<-solve(sqrt(a)%*%(R)%*%sqrt(a))  
    #U_mat=t(Z)%*%a%*%iV%*%t(t((y-mu)))
    y_mu_t=t(t(y-mu)) 
    ZaV<-t(Z)%*%a%*%iV
    ni<-do.call(c,lapply(split(ID$FID,ID$FID),function(x){length(x)}))
    Ub.boot<-matrix(ncol=np,nrow=nrow(U))
    n=length(unique(ID$FID))
    if (pw=="Norm") {
      Ub.boot<-apply(Ub.boot,2, function (x) {
        res.pert<-y_mu_t* rep(rnorm(n), ni) #Normal perturbation
        return(ZaV%*%res.pert)
      })  	
    } else if (pw=="Rade") {
      Ub.boot<-apply(Ub.boot,2, function (x) {
        res.pert<-y_mu_t* rep(sample(c(1,-1),n,replace=T), ni)  #Rademacher distribution
        return(ZaV%*%res.pert)
      })		
    }
    if (Uc==TRUE){Ub.boot<-apply(Ub.boot,1,scale, scale = FALSE)} else if (length(maf)==1) {
      Ub.boot<-t(t(Ub.boot))} else {
        Ub.boot<-t(Ub.boot)} #
    MEANTS=rep(1,11)
    VARTS=rep(1,11)
    DF=rep(1,11)
    for(i in 1:11)
    {
    Ts_boot<-apply(Ub.boot,1,function (x) {t(x)%*%((1-rho[i])*W+rho[i]*tcrossprod(sqrt(W)%*%onev))%*%x })
    
    MEANTS[i]=mean(Ts_boot)
    VARTS[i]=var(Ts_boot)
    DF[i]=12/kurtosis(Ts_boot)
   # var(Ts_boot); var_theory=2*(sum(diag(B%*%t(C)%*%W%*%C%*%B%*%t(C)%*%W%*%C)))
    df=12/kurtosis(Ts_boot)
    if (df<0) {df=100}
    
    mu_Ts<- sum(diag(B%*%t(C)%*%W%*%C)) #trace(BC)  # Bsqrt=mysqrt(B); BC<-Bsqrt%*%t(C)%*%C%*%Bsqrt
    #mu_Ts
    if (df<0) {pval=NA} else {      
      pval[i]=1-pchisq((TS[i]-mean(Ts_boot))*sqrt(2*df)/sqrt(var(Ts_boot))+df, df=df) #pval
    }
    }
  }
  else {
    Ts_boot=mu_Ts=var_theory=df=pval=NA
  }
  pmin=min(pval)  #step 3
  
  Ub.h<-matrix(ncol=NB,nrow=nrow(U))
  Ub.h<-apply(Ub.h,2, function (x) {
    res.pert<-y_mu_t* rep(rnorm(n), ni) #Normal perturbation
    return(ZaV%*%res.pert)
  })  
  pminb=rep(1,NB)
  for(j in 1:NB){
   for(i in 1:11){
  TSh=t(Ub.h[,j])%*%((1-rho[i])*W+rho[i]*tcrossprod(sqrt(W)%*%onev))%*%Ub.h[,j]
  pvalb[i]=1-pchisq((TSh-MEANTS[i])*sqrt(2*DF[i])/sqrt(VARTS[i])+DF[i], df=DF[i])
   }
  pminb[j]=min(pvalb)
  }
  pfinal=sum(pminb<pmin)/NB
 # p1<-results$Qq
 # p2<-pval
 # return(c(p1,p2)) 
 return(pfinal)
}


blockMatrixDiagonal<-function(...){  
  matrixList<-list(...)
  if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
  
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)
  finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
  }
  finalMatrix
} 
##


#square root of a matrix
#mysqrt<-function (V) {
#  #V=B
#  V.eig<-eigen(V)
#  V.eigvalues<-V.eig$values
#  V.eigvalues[-1e-8<V.eigvalues & V.eigvalues<0]=0
#  V.sqrt=V.eig$vectors %*% diag(sqrt(V.eigvalues)) %*% solve(V.eig$vectors)
#  return(V.sqrt)
# }

#square root of a matrix
#from http://realizationsinbiostatistics.blogspot.com/2008/08/matrix-square-roots-in-r_18.html
denman.beavers <- function(mat,maxit=50) {
  stopifnot(nrow(mat) == ncol(mat))
  niter <- 0
  y <- mat
  z <- diag(rep(1,nrow(mat)))
  for (niter in 1:maxit) {
    y.temp <- 0.5*(y+solve(z))
    z <- 0.5*(z+solve(y))
    y <- y.temp
  }
  return(list(sqrt=y,sqrt.inv=z))
}

