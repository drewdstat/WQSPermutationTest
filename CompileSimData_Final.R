#!/usr/local/bin/Rscript

PATH<-"/active/sathyanarayana_s/Lab_Members/DrewDay/WQS_Sim/"

load(paste0(PATH,"WQSPTsim1_FinalPublic_1.RData"))

WQSPTsim1new<-WQSPTsim1$Table
WQSPTsim1new$Rep<-1
WQSPTsim1new$N_null<-0

if(is.null(WQSPTsim1$WQSRHModels)){
  WQSPTsim1new[WQSPTsim1$Table$Model=="WQSBS_RH","N_null"]<-NA
} else {
  WQSPTsim1new[WQSPTsim1$Table$Model=="WQSBS_RH","N_null"]<-
    length(which(sapply(WQSPTsim1$WQSRHModels$gwqslist,is.null)))
}
if(is.null(WQSPTsim1$WQSRSRHModels)){
  WQSPTsim1new[WQSPTsim1$Table$Model=="WQSRS_RH","N_null"]<-NA
} else {
  WQSPTsim1new[WQSPTsim1$Table$Model=="WQSRS_RH","N_null"]<-
    length(which(sapply(WQSPTsim1$WQSRSRHModels$gwqslist,is.null)))
}

WQSPTsim1new[WQSPTsim1$Table$Model=="WQSBS_PT","N_null"]<-
  length(which(is.na(WQSPTsim1$WQSBSPTBetas)))
WQSPTsim1new[WQSPTsim1$Table$Model=="WQSRS_PT","N_null"]<-
  length(which(is.na(WQSPTsim1$WQSRSPTBetas)))

WQSPTsim1$WQSRHModels$wmat<-WQSPTsim1$WQSRHModels$wmat[,paste0("T",1:ncol(WQSPTsim1$WQSRHModels$wmat))]
WQSPTsim1$WQSRSRHModels$wmat<-WQSPTsim1$WQSRSRHModels$wmat[,paste0("T",1:ncol(WQSPTsim1$WQSRSRHModels$wmat))]

AllRHesttable<-cbind(unique(WQSPTsim1new$Seed),"WQSBS_RH",1:nrow(WQSPTsim1$WQSRHModels$coefmat),
                     WQSPTsim1$WQSRHModels$coefmat,WQSPTsim1$WQSRHModels$wmat)
AllRHesttable<-as.data.frame(rbind(AllRHesttable,
                                   cbind(unique(WQSPTsim1new$Seed),"WQSRS_RH",1:nrow(WQSPTsim1$WQSRSRHModels$coefmat),
                                                       WQSPTsim1$WQSRSRHModels$coefmat,WQSPTsim1$WQSRSRHModels$wmat)))
names(AllRHesttable)[1:5]<-c("Seed","Model","Iter","beta0","beta1")

AllPTbetas<-as.data.frame(rbind(cbind(unique(WQSPTsim1new$Seed),"WQSBS_PT",WQSPTsim1$WQSBSPTBetas),
                                cbind(unique(WQSPTsim1new$Seed),"WQSRS_PT",WQSPTsim1$WQSRSPTBetas)))
names(AllPTbetas)<-c("Seed","Model","Permuted_beta1")

for(i in c(2:5)){
  print(i)
  load(paste0(PATH,"WQSPTsim1_FinalPublic_",i,".RData"))
  WQSPTsim1temp<-WQSPTsim1$Table
  WQSPTsim1temp$Rep<-i
  WQSPTsim1temp$N_null<-0
  
  if(!is.null(WQSPTsim1$WQSBSRHModels)&!is.null(WQSPTsim1$WQSRSRHModels)){
    WQSPTsim1temp[WQSPTsim1temp$Table$Model=="WQSBS_RH","N_null"]<-
      length(which(sapply(WQSPTsim1$WQSBSRHModels$gwqslist,is.null)))
    WQSPTsim1temp[WQSPTsim1temp$Table$Model=="WQSRS_RH","N_null"]<-
      length(which(sapply(WQSPTsim1$WQSRSRHModels$gwqslist,is.null)))
    WQSPTsim1$WQSBSRHModels$wmat<-WQSPTsim1$WQSBSRHModels$wmat[,paste0("T",1:ncol(WQSPTsim1$WQSBSRHModels$wmat))]
    WQSPTsim1$WQSRSRHModels$wmat<-WQSPTsim1$WQSRSRHModels$wmat[,paste0("T",1:ncol(WQSPTsim1$WQSRSRHModels$wmat))]
    
    AllRHesttable_temp<-cbind("WQSBS_RH",1:nrow(WQSPTsim1$WQSRHModels$coefmat),
                              WQSPTsim1$WQSRHModels$coefmat,WQSPTsim1$WQSRHModels$wmat)
    AllRHesttable_temp<-as.data.frame(rbind(AllRHesttable_temp,
                                            cbind("WQSRS_RH",1:nrow(WQSPTsim1$WQSRSRHModels$coefmat),
                                                  WQSPTsim1$WQSRSRHModels$coefmat,WQSPTsim1$WQSRSRHModels$wmat)))
    names(AllRHesttable_temp)[1:4]<-c("Model","Iter","beta0","beta1")
    AllRHesttable_temp$Seed<-unique(WQSPTsim1temp$Seed)
    AllRHesttable_temp<-AllRHesttable_temp[,c(ncol(AllRHesttable_temp),1:(ncol(AllRHesttable_temp)-1))]
    AllRHesttable<-rbind(AllRHesttable,AllRHesttable_temp)
  } else if(!is.null(WQSPTsim1$WQSRHModels)&is.null(WQSPTsim1$WQSRSRHModels)){
    WQSPTsim1temp[WQSPTsim1temp$Table$Model=="WQSBS_RH","N_null"]<-
      length(which(sapply(WQSPTsim1$WQSRHModels$gwqslist,is.null)))
    WQSPTsim1temp[WQSPTsim1temp$Table$Model=="WQSRS_RH","N_null"]<-NA
    WQSPTsim1$WQSRHModels$wmat<-WQSPTsim1$WQSRHModels$wmat[,paste0("T",1:ncol(WQSPTsim1$WQSRHModels$wmat))]
    AllRHesttable_temp<-as.data.frame(cbind("WQSBS_RH",1:nrow(WQSPTsim1$WQSRHModels$coefmat),
                                            WQSPTsim1$WQSRHModels$coefmat,WQSPTsim1$WQSRHModels$wmat))
    names(AllRHesttable_temp)[1:4]<-c("Model","Iter","beta0","beta1")
    AllRHesttable_temp$Seed<-unique(WQSPTsim1temp$Seed)
    AllRHesttable_temp<-AllRHesttable_temp[,c(ncol(AllRHesttable_temp),1:(ncol(AllRHesttable_temp)-1))]
    AllRHesttable<-rbind(AllRHesttable,AllRHesttable_temp)
  } else if(is.null(WQSPTsim1$WQSRHModels)&!is.null(WQSPTsim1$WQSRSRHModels)){
    WQSPTsim1temp[WQSPTsim1temp$Table$Model=="WQSRS_RH","N_null"]<-
      length(which(sapply(WQSPTsim1$WQSRSRHModels$gwqslist,is.null)))
    WQSPTsim1temp[WQSPTsim1temp$Table$Model=="WQSBS_RH","N_null"]<-NA
    WQSPTsim1$WQSRSRHModels$wmat<-WQSPTsim1$WQSRSRHModels$wmat[,paste0("T",1:ncol(WQSPTsim1$WQSRSRHModels$wmat))]
    
    AllRHesttable_temp<-as.data.frame(cbind("WQSRS_RH",1:nrow(WQSPTsim1$WQSRSRHModels$coefmat),
                                            WQSPTsim1$WQSRSRHModels$coefmat,WQSPTsim1$WQSRSRHModels$wmat))
    names(AllRHesttable_temp)[1:4]<-c("Model","Iter","beta0","beta1")
    AllRHesttable_temp$Seed<-unique(WQSPTsim1temp$Seed)
    AllRHesttable_temp<-AllRHesttable_temp[,c(ncol(AllRHesttable_temp),1:(ncol(AllRHesttable_temp)-1))]
    AllRHesttable<-rbind(AllRHesttable,AllRHesttable_temp)
  } else {
    WQSPTsim1temp[WQSPTsim1$Table$Model=="WQSBS_RH","N_null"]<-NA
    WQSPTsim1temp[WQSPTsim1$Table$Model=="WQSRS_RH","N_null"]<-NA
  }
  
  if(!is.null(WQSPTsim1$WQSBSPTBetas)){
    WQSPTsim1temp[WQSPTsim1$Table$Model=="WQSBS_PT","N_null"]<-
      length(which(is.na(WQSPTsim1$WQSBSPTBetas)))
    AllPTbetas_temp<-as.data.frame(cbind(unique(WQSPTsim1temp$Seed),"WQSBS_PT",WQSPTsim1$WQSBSPTBetas))
    names(AllPTbetas_temp)<-c("Seed","Model","Permuted_beta1")
    AllPTbetas<-rbind(AllPTbetas,AllPTbetas_temp)
  } else {WQSPTsim1temp[WQSPTsim1$Table$Model=="WQSBS_PT","N_null"]<-NA}
  if(!is.null(WQSPTsim1$WQSRSPTBetas)){
    WQSPTsim1temp[WQSPTsim1$Table$Model=="WQSRS_PT","N_null"]<-
      length(which(is.na(WQSPTsim1$WQSRSPTBetas)))
    AllPTbetas_temp<-as.data.frame(cbind(unique(WQSPTsim1temp$Seed),"WQSRS_PT",WQSPTsim1$WQSBSPTBetas))
    names(AllPTbetas_temp)<-c("Seed","Model","Permuted_beta1")
    AllPTbetas<-rbind(AllPTbetas,AllPTbetas_temp)
  } else {WQSPTsim1temp[WQSPTsim1$Table$Model=="WQSRS_PT","N_null"]<-NA}
  
  WQSPTsim1new<-rbind(WQSPTsim1new,WQSPTsim1temp)
}

WQSPTsim1new$Parameter2<-ifelse(WQSPTsim1new$Class=="Beta"&WQSPTsim1new$Parameter%in%
                                  c("beta0","beta1","WQSbeta","bpos","bneg"),
                                as.character(WQSPTsim1new$Parameter),
                                ifelse(grepl("C|gamma",WQSPTsim1new$Parameter),
                                       "CovBeta","Weight"))
as.numeric.factor<-function(x) as.numeric(levels(x))[x]
AllPTbetas$Permuted_beta1<-as.numeric.factor(AllPTbetas$Permuted_beta1)
AllRHesttable[,3:ncol(AllRHesttable)]<-lapply(AllRHesttable[,3:ncol(AllRHesttable)],as.numeric.factor)
Allmodssim1_finaltest<-list(Main=WQSPTsim1new,AllPTbetas=AllPTbetas,
                                         AllRHest=AllRHesttable)
save(Allmodssim1_finaltest,
     file=paste0(PATH,"Allmodssim1_finaltest.RData"))
