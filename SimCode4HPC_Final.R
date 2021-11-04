#!/usr/local/bin/Rscript

PATH<-"/MyFolder/" #Put your own file path here

library(mvtnorm)
library(extraDistr)
library(ggplot2)
library(pbapply)
library(plyr)
library(dplyr)
library(gWQS)

args<-commandArgs(trailingOnly=TRUE)
# ---------------------------------------------------------------------------------------------------------

source(paste0(PATH,"permtest_v304_final.R"))
source(paste0(PATH,"WQSsim.R"))

mytruegammas<-c(-0.63,0.18,-0.84,1.60,0.33,rep(0,5))

Sim_JustPT<-function(simwide_b1pos=T,nmixs=10,ptruewts=1,nobss=500,ncovrts=10,ptruecovrts=0.5,
                     corrstructs=0,epss=1,truewqsbeta=0.3,truebeta0=2,
                     truewts=c(rep(0.15,5),rep(0.05,5)),truegamma=mytruegammas,
                     boots=100,nrs=100,quantile=T,nq=5,pt.nrep=200,myplan="multicore",rephold.wqs=T,
                     rephold.wqsrs=T,rh.nrep=100,wqsPT=T,wqsrsPT=T,includeqgc=T,seed=NULL){
  if(is.null(truewqsbeta)==T){
    mytruewqsbeta=NULL
  } else {
    mytruewqsbeta=truewqsbeta
  }
  if(is.null(truebeta0)==T){
    mytruebeta0=NULL
  } else {
    mytruebeta0=truebeta0
  }
  if(is.null(truewts)==T){
    mytruewts=NULL
  } else {
    mytruewts=truewts
  }
  if(is.null(truegamma)==T){
    mytruegamma=NULL
  } else {
    mytruegamma=truegamma
  }
  if(is.null(seed)){
    myseed<-sample(1:1E9,1)
  } else {
    myseed<-seed
  }
  
  newsim<-WQSsim(nmix=nmixs,ncovrt=ncovrts,nobs=nobss,ntruewts=(nmixs*ptruewts),ntruecovrt=(ncovrts*ptruecovrts),
                 corrstruct=corrstructs,eps=epss,truewqsbeta=mytruewqsbeta,truebeta0=mytruebeta0,
                 truewts=mytruewts,truegamma=mytruegamma,constrdir="none",seed=myseed,q=nq)
  
  if(ncovrts!=0){
    form1<-formula(paste0("y~wqs+",paste(paste0("C",1:ncovrts),collapse="+")))
  } else {
    form1<-formula(paste0("y~wqs"))
  }
  
  mix_name<-names(newsim$Data)[grep("T",names(newsim$Data))]
  
  if(includeqgc==T){
    library(qgcomp)
    if(ncovrts!=0){
      qgcform1<-formula(paste0("y~",paste(paste0("T",1:nmixs),collapse="+"),"+",
                               paste(paste0("C",1:ncovrts),collapse="+")))
    } else {
      qgcform1<-formula(paste0("y~."))
    }
    qgc.nb<-qgcomp.noboot(qgcform1,data=newsim$Data,
                          expnms=names(newsim$Data)[grep("T",names(newsim$Data))],
                          q=nq)
    qgc.b<-qgcomp.boot(qgcform1,data=newsim$Data,
                       expnms=names(newsim$Data)[grep("T",names(newsim$Data))],
                       q=nq,B=boots*10)
    print("qgcs complete")
    extract_qgc<-function(qgcmod,simdat=newsim,qgcmodname="QGC_Boot"){
      se_comb<-function (expnms, covmat, grad = NULL){
        if (!is.matrix(covmat)) {
          nm <- names(covmat)
          covmat = matrix(covmat)
          colnames(covmat) <- nm
        }
        weightvec <- rep(0, dim(covmat)[1])
        if (is.null(grad)) 
          weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- 1
        if (!is.null(grad)) 
          weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- grad
        var <- weightvec %*% covmat %*% weightvec
        sqrt(var)[1, , drop = TRUE]
      }
      outmat<-as.data.frame(matrix(NA,(length(simdat$coef)+nrow(simdat$weights)+2),7))
      names(outmat)<-c("Model","Parameter","True_value","mean","CI_2.5","CI_97.5","p-value")
      outmat$Model<-qgcmodname
      outmat$Parameter<-c(names(simdat$coef),as.character(simdat$weights$mix_name),"bpos","bneg")
      outmat$True_value<-c(simdat$coef,simdat$coef[2]*simdat$weights$true_weight,simdat$coef[2],0)
      
      if("msmfit"%in%names(qgcmod)){ #this is for qgcomp.boot
        psum<-sum(qgcmod$fit$coef[which(dimnames(summary(qgcmod$fit)$coef)[[1]]%in%qgcmod$expnms&
                                          qgcmod$fit$coef>0)])
        nsum<-sum(qgcmod$fit$coef[which(dimnames(summary(qgcmod$fit)$coef)[[1]]%in%qgcmod$expnms&
                                          qgcmod$fit$coef<0)])
        outmat$mean<-c(qgcmod$fit$coef[1],qgcmod$psi,
                       qgcmod$fit$coef[grep("C",names(qgcmod$fit$coef))],
                       qgcmod$fit$coef[grep("T",names(qgcmod$fit$coef))],psum,nsum)
        outmat[,c("CI_2.5","CI_97.5")]<-
          rbind(confint(qgcmod$fit)[1,],qgcmod$ci,
                confint(qgcmod$fit)[grep("C",dimnames(confint(qgcmod$fit))[[1]]),],
                confint(qgcmod$fit)[grep("T",dimnames(confint(qgcmod$fit))[[1]]),],
                c(NA,NA),c(NA,NA))
        outmat[,"p-value"]<-
          c(summary(qgcmod$fit)$coef[1,4],qgcmod$pval[length(qgcmod$pval)],
            summary(qgcmod$fit)$coef[grep("C",dimnames(summary(qgcmod$fit)$coef)[[1]]),4],
            summary(qgcmod$fit)$coef[grep("T",dimnames(summary(qgcmod$fit)$coef)[[1]]),4],NA,NA)
      } else { #this is for qgcomp.noboot
        outmat$mean<-c(qgcmod$fit$coef[1],qgcmod$psi,
                       qgcmod$fit$coef[grep("C",names(qgcmod$fit$coef))],
                       qgcmod$fit$coef[grep("T",names(qgcmod$fit$coef))],qgcmod$pos.psi,qgcmod$neg.psi)
        outmat[,c("CI_2.5","CI_97.5")]<-
          rbind(confint(qgcmod$fit)[1,],qgcmod$ci,
                confint(qgcmod$fit)[grep("C",dimnames(confint(qgcmod$fit))[[1]]),],
                confint(qgcmod$fit)[grep("T",dimnames(confint(qgcmod$fit))[[1]]),],
                c(NA,NA),c(NA,NA))
        outmat$'p-value'<-c(summary(qgcmod$fit)$coef[1,4],qgcmod$pval[length(qgcmod$pval)],
                            summary(qgcmod$fit)$coef[grep("C",dimnames(summary(qgcmod$fit)$coef)[[1]]),4],
                            summary(qgcmod$fit)$coef[grep("T",dimnames(summary(qgcmod$fit)$coef)[[1]]),4],NA,NA)
      }
      outmat$Class<-c(rep("Beta",length(simdat$betas)),rep("Weight",nrow(simdat$weights)),rep("Beta",2))
      outmat$Capture<-
        ifelse(outmat$True_value>outmat$CI_2.5&outmat$True_value<outmat$CI_97.5,1,0)
      outmat$Error<-outmat$mean-outmat$True_value
      outmat$False_positive<-ifelse(outmat$True_value==0&
                                      outmat$CI_2.5<0&outmat$CI_97.5<0,1,
                                    ifelse(outmat$True_value==0&outmat$CI_2.5>0&outmat$CI_97.5>0,1,
                                           ifelse(is.na(outmat$CI_2.5),NA,0)))
      outmat$False_negative<-ifelse(outmat$True_value!=0&
                                      outmat$CI_2.5<0&outmat$CI_97.5>0,1,
                                    ifelse(is.na(outmat$CI_2.5),NA,0))
      return(outmat)
    }
    qgcextract1<-extract_qgc(qgc.nb,simdat=newsim,qgcmodname="QGC_Noboot")
    qgcextract2<-extract_qgc(qgc.b,simdat=newsim,qgcmodname="QGC_Boot")
    qgcextracts<-rbind(qgcextract1,qgcextract2)
  }
  
  newwqs<-tryCatch(
    gwqs(formula=form1,mix_name=names(newsim$Data)[grep("T",names(newsim$Data))],data=newsim$Data,na.action=na.exclude,q=nq,validation=0,b=boots,plan_strategy=myplan,b1_pos=simwide_b1pos),
    error=function(e) NULL)
  newwqs2<-tryCatch(
    gwqs(formula=form1,mix_name=names(newsim$Data)[grep("T",names(newsim$Data))],data=newsim$Data,na.action=na.exclude,q=nq,validation=0.6,b=boots,plan_strategy=myplan,b1_pos=simwide_b1pos),
    error=function(e) NULL)
  newwqs3<-tryCatch(
    gwqs(formula=form1,mix_name=names(newsim$Data)[grep("T",names(newsim$Data))],data=newsim$Data,na.action=na.exclude,q=nq,validation=0,b=nrs,rs=T,plan_strategy=myplan,b1_pos=simwide_b1pos),
    error=function(e) NULL)
  newwqs4<-tryCatch(
    gwqs(formula=form1,mix_name=names(newsim$Data)[grep("T",names(newsim$Data))],data=newsim$Data,na.action=na.exclude,q=nq,validation=0.6,b=nrs,rs=T,plan_strategy=myplan,b1_pos=simwide_b1pos),
    error=function(e) NULL)
  #save.image(file="IntermediateSimCode4HPC_FinalCheck.RData") pick up here
  WQSoutput<-function(model,sim,pt=T,modtype="WQS_Nosplit",pt.nr,bts,RS=F,bplus=T){
    if(is.null(model)|is.null(model$final_weights)){
      outmat<-as.data.frame(matrix(NA,(length(sim$coef)+nrow(sim$weights)),12))
      names(outmat)<-c("Model","Parameter","True_value","mean","CI_2.5","CI_97.5","p-value",
                       "Class","Capture","Error","False_positive","False_negative")
      outmat$Model<-modtype
      outmat$Parameter<-c(names(sim$coef),as.character(sim$weights$mix_name))
      outmat$True_value<-c(sim$coef,sim$weights$true_weight)
      outmat$Class<-c(rep("Beta",length(sim$coef)),rep("Weight",nrow(sim$weights)))
      #outmat$False_positive<-0
      #outmat$False_negative<-ifelse(outmat$True_value==0,0,1)
      #outmat$Error<-0-outmat$True_value
      if(pt){
        outmat2<-outmat
        if(RS){outmat2$Model<-"WQSRS_PT"} else {outmat2$Model<-"WQS_PT"}
        outmat<-rbind(outmat,outmat2)
        betamat<-NULL
      }
    } else {
      outmat<-as.data.frame(matrix(NA,(length(sim$coef)+nrow(sim$weights)),7))
      names(outmat)<-c("Model","Parameter","True_value","mean","CI_2.5","CI_97.5","p-value")
      outmat$Model<-modtype
      outmat$Parameter<-c(names(sim$coef),as.character(sim$weights$mix_name))
      outmat$True_value<-c(sim$coef,sim$weights$true_weight)
      if(!"mix_name"%in%names(model$final_weights)){ 
        names(model$final_weights)<-c("mix_name","mean_weight")
      }
      wtmat<-merge(sim$weights,model$final_weights,by="mix_name",all=T)
      wtmat$mix_name<-factor(wtmat$mix_name,levels=as.character(sim$weights$mix_name))
      wtmat<-wtmat[order(wtmat$mix_name),]
      
      outmat$mean<-c(model$fit$coef,wtmat$mean_weight)
      outmat[1:length(sim$coef),c("CI_2.5","CI_97.5")]<-confint(model$fit)
      outmat[1:length(sim$coef),"p-value"]<-summary(model$fit)$coef[,4]
      outmat$Class<-c(rep("Beta",length(sim$coef)),rep("Weight",nrow(sim$weights)))
      
      if(pt==T){
        if(RS==T){
          newPT<-perm.test(model,niter=pt.nr,boots=bts,b1pos=bplus,RS=T)
          outmat2<-outmat
          outmat2$Model<-"WQSRS_PT"
          outmat2[2,"p-value"]<-newPT$pval
          outmat<-rbind(outmat,outmat2)
          betamat<-newPT$betas
        } else {
          newPT<-perm.test(model,niter=pt.nr,b1pos=bplus,boots=bts)
          outmat2<-outmat
          outmat2$Model<-"WQS_PT"
          outmat2[2,"p-value"]<-newPT$pval
          outmat<-rbind(outmat,outmat2)
          betamat<-newPT$betas
        }
      } 
      
      outmat$Capture<-
        ifelse(outmat$True_value>outmat$CI_2.5&outmat$True_value<outmat$CI_97.5,1,0)
      outmat$Error<-outmat$mean-outmat$True_value
      outmat$False_positive<-ifelse(outmat$Class=="Beta"&outmat$True_value==0&outmat$'p-value'<0.05,1,
                                    ifelse(outmat$Class=="Weight",NA,0))
      outmat$False_negative<-ifelse(outmat$Class=="Beta"&outmat$True_value!=0&outmat$'p-value'>=0.05,1,
                                    ifelse(outmat$Class=="Weight",NA,0))
    }
    if(pt){
      retlist<-list(Table=outmat,Betas=betamat)
    } else {
      retlist<-outmat
    }
    return(retlist)
  }
  
  if(wqsPT==T){
    wout1<-WQSoutput(model=newwqs,sim=newsim,pt=T,modtype="WQS_Nosplit",pt.nr=pt.nrep,bts=boots,bplus=simwide_b1pos)
    myout<-wout1$Table
    WQSPTBetas<-wout1$Betas
  } else {
    myout<-WQSoutput(model=newwqs,sim=newsim,pt=F,modtype="WQS_Nosplit",pt.nr=pt.nrep,bts=boots,bplus=simwide_b1pos)
  }
  myout<-rbind(myout,WQSoutput(model=newwqs2,sim=newsim,pt=F,modtype="WQS_Split",pt.nr=0,bts=boots,bplus=simwide_b1pos))
  if(wqsrsPT==T){
    wout2<-WQSoutput(model=newwqs3,sim=newsim,pt=T,modtype="WQSRS_Nosplit",pt.nr=pt.nrep,
                     bts=nrs,RS=T,bplus=simwide_b1pos)
    myout<-rbind(myout,wout2$Table)
    WQSRSPTBetas<-wout2$Betas
  } else {
    myout<-rbind(myout,WQSoutput(model=newwqs3,sim=newsim,pt=F,modtype="WQSRS_Nosplit",pt.nr=pt.nrep,
                                 bts=nrs,RS=T,bplus=simwide_b1pos))
  }
  myout<-rbind(myout,WQSoutput(model=newwqs4,sim=newsim,pt=F,modtype="WQSRS_Split",pt.nr=0,
                               bts=nrs,RS=T,bplus=simwide_b1pos))
  
  processrh<-function(rhmodel,modelname="WQS_RH"){
    if(is.null(rhmodel)){
      outmat<-data.frame(Model=modelname,
                         Parameter=c("beta0","beta1",paste0("gamma",1:ncovrts),paste0("T",1:nmixs)),
                         True_value=c(truebeta0,truewqsbeta,truegamma,truewts),mean=NA,CI_2.5=NA,CI_97.5=NA)
      outmat[,"p-value"]<-NA
      outmat$Class<-ifelse(grepl("beta|gamma",outmat$Parameter),"Beta","Weight")
      outmat[,c("Capture","Error","False_positive","False_negative")]<-NA
    } else {
      outmat<-data.frame(Model=modelname,
                         Parameter=c("beta0","beta1",paste0("gamma",1:ncovrts),paste0("T",1:nmixs)),
                         True_value=c(truebeta0,truewqsbeta,truegamma,truewts))
      outmat[,c("mean","CI_2.5","CI_97.5")]<-
        rbind(as.matrix(rhmodel$fit$coefficients[,c(1,3:4)]),
              as.matrix(rhmodel$final_weights[match(paste0("T",1:nmixs),rhmodel$final_weights$mix_name),
                                              c(2:4)]))
      outmat$'p-value'<-NA
      outmat$Class<-ifelse(grepl("beta|gamma",outmat$Parameter),"Beta","Weight")
      outmat$Capture<-
        ifelse(outmat$True_value>outmat$CI_2.5&outmat$True_value<outmat$CI_97.5,1,0)
      outmat$Error<-outmat$mean-outmat$True_value
      outmat$False_positive<-ifelse(outmat$Class=="Beta"&outmat$True_value!=0,0,
                                    ifelse(outmat$Class=="Beta"&outmat$True_value==0&outmat$CI_97.5>0&
                                             outmat$CI_2.5<0,0,ifelse(outmat$Class=="Weight",NA,1)))
      outmat$False_negative<-ifelse(outmat$Class=="Beta"&outmat$True_value!=0&outmat$CI_97.5>0&
                                      outmat$CI_2.5<0,1,ifelse(outmat$Class=="Weight",NA,0))
    }
    return(outmat)
  }
  
  if(rephold.wqs==T&rephold.wqsrs==T){
    rep.hold<-tryCatch({
      gwqsrh(formula=formula(paste0("y~wqs+",paste(paste0("C",1:10),collapse="+"))),
               mix_name=names(newsim$Data)[grep("T",names(newsim$Data))],
               data=newsim$Data,na.action=na.exclude,q=nq,validation=0.6,b=boots,rs=F,
               plan_strategy=myplan,b1_pos=simwide_b1pos,rh=rh.nrep)
    },error=function(e) NULL)
    rep.hold2<-tryCatch({
      gwqsrh(formula=formula(paste0("y~wqs+",paste(paste0("C",1:10),collapse="+"))),
               mix_name=names(newsim$Data)[grep("T",names(newsim$Data))],
               data=newsim$Data,na.action=na.exclude,q=nq,validation=0.6,b=nrs,rs=T,
               plan_strategy=myplan,b1_pos=simwide_b1pos,rh=rh.nrep)
    },error=function(e) NULL)
    myout<-rbind(myout,processrh(rep.hold))
    myout<-rbind(myout,processrh(rep.hold2,modelname="WQSRS_RH"))
  } else if(rephold.wqs==T&rephold.wqsrs==F){
    rep.hold<-tryCatch({
      gwqsrh(formula=form1,mix_name=names(newsim$Data)[grep("T",names(newsim$Data))],
               data=newsim$Data,na.action=na.exclude,q=nq,validation=0.6,b=boots,rs=F,
               plan_strategy=myplan,b1_pos=simwide_b1pos,rh=rh.nrep)
    },error=function(e) NULL)
    myout<-rbind(myout,processrh(rep.hold))
  } else if(rephold.wqs==F&rephold.wqsrs==T){
    rep.hold2<-tryCatch({
      gwqsrh(formula=form1,mix_name=names(newsim$Data)[grep("T",names(newsim$Data))],
               data=newsim$Data,na.action=na.exclude,q=nq,validation=0.6,b=nrs,rs=T,
               plan_strategy=myplan,b1_pos=simwide_b1pos,rh=rh.nrep)
    },error=function(e) NULL)
    myout<-rbind(myout,processrh(rep.hold2,modelname="WQSRS_RH"))
  }
  
  if(includeqgc==T){
    myout<-rbind(myout,qgcextracts)
  }
  
  myout$Rep<-args[1]
  myout$Seed<-myseed
  myout$Nmix<-nmixs
  myout$Ntruewt<-ptruewts*nmixs
  myout$Nobs<-nobss
  myout$Ncovrt<-ncovrts
  myout$Ntruecovrt<-ptruecovrts*ncovrts
  if(is.matrix(corrstructs)){
    myout$Corrstruct<-"Real VC"
  } else {
    myout$Corrstruct<-corrstructs
  }
  myout$Quantile<-quantile
  myout$Eps<-epss
  myout$True_beta1<-mytruewqsbeta
  myout$True_beta0<-mytruebeta0
  
  outlist<-list(Table=myout)
  if(wqsPT){
    outlist[["WQSPTBetas"]]<-WQSPTBetas
  } else {
    outlist[["WQSPTBetas"]]<-NULL
  }
  if(wqsrsPT){
    outlist[["WQSRSPTBetas"]]<-WQSRSPTBetas
  } else {
    outlist[["WQSRSPTBetas"]]<-NULL
  }
  if(rephold.wqs){
    outlist[["WQSRHModels"]]<-rep.hold
  } else {
    outlist[["WQSRHModels"]]<-NULL
  }
  if(rephold.wqsrs){
    outlist[["WQSRSRHModels"]]<-rep.hold2
  } else {
    outlist[["WQSRSRHModels"]]<-NULL
  }
  return(outlist)
}

#True beta1 = 0.3, Predictor Correlation = Uncorrelated
WQSPTsim1<-Sim_JustPT(truewqsbeta=0.3,corrstructs=0)
save(WQSPTsim1,file=paste0(PATH,"WQSPTsim1_UncorrNZ_",args[1],".RData"))

#True beta1 = 0, Predictor Correlation = Uncorrelated
WQSPTsim1<-Sim_JustPT(truewqsbeta=0,corrstructs=0)
save(WQSPTsim1,file=paste0(PATH,"WQSPTsim1_UncorrNZ_",args[1],".RData"))

#True beta1 = 0.2, Predictor Correlation = Real Variance-Covariance Matrix
load(paste0(PATH,"examplephxvcov.RData")) #derived from TIDES cohort data (see Methods)
WQSPTsim1<-Sim_JustPT(truewqsbeta=0.2,corrstructs=cov2cor(examplephxvcov))
save(WQSPTsim1,file=paste0(PATH,"WQSPTsim1_UncorrNZ_",args[1],".RData"))

#True beta1 = 0, Predictor Correlation = Real Variance-Covariance Matrix
WQSPTsim1<-Sim_JustPT(truewqsbeta=0,corrstructs=cov2cor(examplephxvcov))
save(WQSPTsim1,file=paste0(PATH,"WQSPTsim1_CorrZ_",args[1],".RData"))

