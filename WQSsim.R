WQSsim<-function(nmix=10,ncovrt=10,nobs=500,ntruewts=10,ntruecovrt=5,corrstruct=0,eps=1,
                 truewqsbeta=NULL,truebeta0=NULL,truewts=NULL,truegamma=NULL,
                 constrdir="none",seed=101,q=10,binix=F,binixwt=F,ixbeta=NULL,bindiff=NULL){
  library(mvtnorm)
  library(extraDistr)
  
  if(length(corrstruct)==1){
    Rho<-diag(nmix+ncovrt)
    Rho[upper.tri(Rho)]<-Rho[lower.tri(Rho)]<-corrstruct
  } else {Rho<-corrstruct}
  
  if(binixwt==F){
    weights<-rep(0,nmix)
    if(is.null(truewts)){
      set.seed(seed)
      truewts<-extraDistr::rdirichlet(1,rep(1,ntruewts))
      weights[1:ntruewts]<-truewts
    } else {
      if(length(truewts)==nmix&sum(abs(truewts))!=1){
        truewts<-truewts/sum(truewts)
      } 
      if(length(truewts)<nmix){
        weights[1:length(truewts)]<-truewts
        weights[(length(truewts)+1):nmix]<-(1-sum(truewts))/(nmix-length(truewts))
      } else {
        weights[1:length(truewts)]<-truewts
      }
    }
  } else {
    weights<-rep(0,nmix*2)
    if(is.null(truewts)){
      set.seed(seed)
      truewts<-extraDistr::rdirichlet(1,rep(1,ntruewts))
      weights[1:ntruewts]<-truewts
    } else {
      if(length(truewts)==(nmix*2)&sum(abs(truewts))!=1){
        truewts<-truewts/sum(truewts)
      } 
      if(length(truewts)!=(nmix*2)){
        weights[1:length(truewts)]<-truewts
        weights[(length(truewts)+1):(nmix*2)]<-(1-sum(truewts))/((nmix*2)-length(truewts))
      } else {
        weights[1:length(truewts)]<-truewts
      }
    }
  }
  
  if(round(sum(weights),3)!=1.0){
    warning(print(paste0("weights add up to ",sum(weights))))
  }
  
  set.seed(seed)
  Xmat<-rmvnorm(nobs,mean=rep(0,nmix+ncovrt),sigma=Rho)
  if(is.null(q)){
    Xmatquant<-Xmat
  } else {
    Xmatquant<-Xmat
    Xmatquant[,1:nmix]<-apply(Xmatquant[,1:nmix],2,FUN=function(x){
      as.numeric(as.character(cut(x,breaks=quantile(x,probs=seq(0,1,by=(1/q))),
                                  include.lowest = T,labels=0:(q-1))))
    })
  }
  
  if(binix==T){
    if(is.null(bindiff)) {
      set.seed(seed)
      bindiff<-rnorm(1)
    }
    set.seed(seed)
    Sex<-sample(c(0,1),nobs,replace=T)
    Xmat<-cbind(Xmat,Sex)
    dimnames(Xmat)[[2]][ncol(Xmat)]<-"Sex"
    Xmatquant<-cbind(Xmatquant,Sex)
    dimnames(Xmatquant)[[2]][ncol(Xmatquant)]<-"Sex"
  }
  
  if(ncovrt<ntruecovrt){ntruecovrt<-ncovrt}
  if(!is.null(truegamma)){
    if(length(truegamma)==1){
      covrtbetas<-rep(truegamma,ncovrt)
    } else {
      covrtbetas<-truegamma
    }
  } else {
    set.seed(seed)
    covrtbetas<-c(rnorm(ntruecovrt),rep(0,length=ncovrt-ntruecovrt))
  }
  
  set.seed(seed)
  if(!is.null(truebeta0)){
    beta0<-truebeta0
  } else {
    beta0<-rnorm(1)
  }

  if(!is.null(truewqsbeta)){
    wqsbeta<-truewqsbeta
  } else {
    set.seed(seed)
    if(constrdir=="positive"){
      wqsbeta<-extraDistr::rhnorm(1)
    } else if(constrdir=="negative"){
      wqsbeta<-extraDistr::rhnorm(1)*-1
    } else {
      wqsbeta<-rnorm(1)
    }
  }
  
  if(binix==F){
    wqs<-Xmatquant[,1:nmix]%*%weights
    if(ncovrt>0){
      modmat<-cbind(1,wqs,Xmat[,c((nmix+1):(nmix+ncovrt))])
      dimnames(modmat)[[2]]<-c("Intercept","wqs",paste0("C",1:ncovrt))
      betas<-c(beta0,wqsbeta,covrtbetas)
      names(betas)<-c("beta0","beta1",paste0("gamma",1:ncovrt))
    } else {
      modmat<-cbind(1,wqs)
      dimnames(modmat)[[2]]<-dimnames(modmatq)[[2]]<-c("Intercept","wqs")
      betas<-c(beta0,wqsbeta)
      names(betas)<-c("beta0","beta1")
    }
  } else {
    if(is.null(ixbeta)){
      set.seed(seed)
      ixbeta<-rnorm(1)
    }
    if(binixwt==T){
      sex0rows<-which(Sex==0)
      sex1rows<-which(Sex==1)
      wqs<-rep(0,nobs)
      wqs[sex0rows]<-Xmatquant[sex0rows,1:nmix]%*%weights[1:nmix]
      wqs[sex1rows]<-Xmatquant[sex1rows,1:nmix]%*%weights[(nmix+1):(nmix*2)]
    } else {
      wqs<-Xmatquant[,1:nmix]%*%weights
    }
    
    if(ncovrt>0){
      modmat<-cbind(1,wqs,Sex,Xmatquant[,c((nmix+1):(nmix+ncovrt))])
      dimnames(modmat)[[2]]<-c("Intercept","wqs","Sex",paste0("C",1:ncovrt))
      modmat<-cbind(modmat,modmat[,"wqs"]*modmat[,"Sex"])
      dimnames(modmat)[[2]][ncol(modmat)]<-"wqs_Sex"
      
      betas<-c(beta0,wqsbeta,bindiff,covrtbetas,ixbeta)
      names(betas)<-c("beta0","beta1","sexbeta",paste0("gamma",1:ncovrt),"ixbeta")
    } else {
      modmat<-cbind(1,wqs,Sex)
      dimnames(modmat)[[2]]<-c("Intercept","wqs","Sex")
      modmat<-cbind(modmat,modmat[,"wqs"]*modmat[,"Sex"])
      dimnames(modmat)[[2]][ncol(modmat)]<-"wqs_Sex"
      
      betas<-c(beta0,wqsbeta,bindiff,ixbeta)
      names(betas)<-c("beta0","beta1","sexbeta","ixbeta")
    }
  }
  yhat<-modmat%*%betas
  set.seed(seed)
  epsilon<-rnorm(nobs,sd=eps)
  y<-yhat+epsilon
  Data<-data.frame(cbind(y,Xmat))
  
  if(binix==T){
    if(ncovrt>0){
      names(Data)<-c("y",paste0("T",1:nmix),paste0("C",1:ncovrt),"Sex")
      colnames(Xmatquant)<-c(paste0("T",1:nmix),paste0("C",1:ncovrt),"Sex")
      colnames(Xmat)<-c(paste0("T",1:nmix),paste0("C",1:ncovrt),"Sex")
    } else {
      names(Data)<-c("y",paste0("T",1:nmix),"Sex")
      colnames(Xmatquant)<-c(paste0("T",1:nmix),"Sex")
      colnames(Xmat)<-c(paste0("T",1:nmix),"Sex")
    }
  } else {
    if(ncovrt>0){
      names(Data)<-c("y",paste0("T",1:nmix),paste0("C",1:ncovrt))
      colnames(Xmatquant)<-c(paste0("T",1:nmix),paste0("C",1:ncovrt))
      colnames(Xmat)<-c(paste0("T",1:nmix),paste0("C",1:ncovrt))
    } else {
      names(Data)<-c("y",paste0("T",1:nmix))
      colnames(Xmatquant)<-c(paste0("T",1:nmix))
      colnames(Xmat)<-c(paste0("T",1:nmix))
    }
  }
  if(binixwt==T){
    wtmat<-data.frame(mix_name=rep(paste0("T",1:nmix),2),true_weight=weights,Sex=c(rep(0,nmix),rep(1,nmix)))
  } else {
    wtmat<-data.frame(mix_name=paste0("T",1:nmix),true_weight=weights)
  }
  retlist<-list(weights=wtmat,coef=betas,Data=Data,yhat=yhat,wqs=wqs,modmat=modmat,
                Xq=Xmatquant)
  return(retlist)
}
