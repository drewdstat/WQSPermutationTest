perm.test<-function(model,niter=200,boots=200,b1pos=T,RS=F,myplan="multicore",
                    returnbetas=T){
  library(pbapply)
  
  if(class(model)=="gwqs"){
    mm<-model$fit
    Data<-model$data[model$vindex,-which(names(model$data)%in%c("wqs","wghts"))]
    modclass<-"gwqs"
    if(!is.null(model$qi)){
      nq<-max(sapply(model$qi,length))-1
    } else {
      nq<-NULL 
    }
  } else stop("'model' must be of class 'gwqs' (see gWQS package).")
  
  formchar<-as.character(formula(mm))
  if(!is.null(model$stratified)|grepl("wqs:",formchar[3],fixed=T)){
    stop("This permutation test is not yet set up to accomodate stratified weights or 
         WQS interaction terms.")
  } 
  yname<-as.character(formula(mm))[2]
  if(length(mm$coef)>2){
    fit.partial<-lm(formula(paste0(formchar[2],formchar[1],gsub("wqs + ","",formchar[3],fixed=T))),
                    data=Data)
    partial.yhat<-predict(fit.partial)
    partial.resid<-resid(fit.partial)
    reorgmat<-matrix(NA,dim(Data)[1],niter)
    reorgmat<-apply(reorgmat,2,function(x) partial.yhat+sample(partial.resid,replace=F))
  } else {
    reorgmat<-matrix(NA,dim(Data)[1],niter)
    reorgmat<-apply(reorgmat,2,function(x) sample(Data[,yname])) 
  }
  getbetas<-function(x){ 
    mix_name<-names(model$bres)[names(model$bres)%in%model$final_weights$mix_name]
    newDat<-Data
    newDat[,yname]<-x
    names(newDat)<-c(names(Data))
    formchar<-as.character(formula(mm))
    if(length(mm$coef)>2){
      form1<-formula(paste0(formchar[2],formchar[1],formchar[3]))
    } else {
      form1<-formula(paste0(formchar[2],formchar[1],"wqs"))
    }
    if(RS==T){
      gwqs1<-tryCatch({
        suppressWarnings(gwqs(formula=form1,data=newDat,
                              mix_name=names(model$bres)[names(model$bres)%in%
                                                           model$final_weights$mix_name],q=nq,
                              b=boots,rs=T,validation=0,plan_strategy=myplan,b1_pos=b1pos))
      },error=function(e) NULL, warning=function(e) message("WQSRS failed"))
    } else {
      gwqs1<-tryCatch({
        suppressWarnings(gwqs(formula=form1,data=newDat,
                              mix_name=names(model$bres)[names(model$bres)%in%
                                                           model$final_weights$mix_name],q=nq,
                              b=boots,validation=0,plan_strategy=myplan,b1_pos=b1pos))
      },error=function(e) NULL, warning=function(e) message("WQS failed"))
    }
    if(is.null(gwqs1)) lm1<-NULL else lm1<-gwqs1$fit
    if(is.null(lm1)){retvec<-NA} else {retvec<-lm1$coef[2]} 
    return(retvec)
  }
  
  pboptions(type="timer")
  betas<-pbapply(reorgmat,2,getbetas)
  if(any(is.na(betas))){
    print(paste0(length(which(is.na(betas)))," failed model attempts"))}
  pval<-function(x,true,posb1=b1pos){
    if(posb1){
      length(which(x>true))
    } else {
      length(which(x<true))
    }
  }
  if(returnbetas){
    retlist<-list(pval=pval(betas,mm$coef[2])/length(betas),testbeta1=mm$coef[2],betas=betas)
  } else {
    retlist<-list(pval=pval(betas,mm$coef[2])/length(betas),testbeta1=mm$coef[2])
  }
  return(retlist) 
}
