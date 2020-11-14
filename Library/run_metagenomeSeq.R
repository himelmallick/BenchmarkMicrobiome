#################
# metagenomeSeq #
#################

###########################
# Load Essential Packages #
###########################

# pacman, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
load_essential_packages()

###########################################
# Load Dedicated Method-specific Packages #
###########################################

pacman::p_load('reshape2')

if(! require("metagenomeSeq")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("metagenomeSeq")
}
suppressPackageStartupMessages(library(metagenomeSeq))

##################################
# Fit metagenomeSeq To A Dataset #
##################################

fit.metagenomeSeq <- function(features, 
                              metadata, 
                              libSize, 
                              ID, 
                              transformation,
                              multiple_qvalues){
  
  #########################
  # Transformation if any #
  #########################
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default metagenomeSeq model. Use NONE.')
  
  ###################################
  # Standard metagenomeSeq pipeline #
  ###################################
  
  design <- model.matrix(~., data=metadata)
  count_table <- t(features) 
  mgsdata <- newMRexperiment(counts = count_table)
  mgsp <- cumNormStat(mgsdata)
  mgsdata <- cumNorm(mgsdata, mgsp)
  
  #####################
  # Fit metagenomeSeq # 
  #####################
  
  ############################
  # Random Effect Adjustment #
  ############################
  
  if(!length(ID)==length(unique(ID))){
    fit <- fitZig(obj=mgsdata,mod=design, useMixedModel=TRUE,block=ID)
  } else{
    fit <- fitZig(obj=mgsdata,mod=design)
  }
  
  ###################
  # Combine results #
  ###################
  
  if(dim(metadata)[2]>1){
    coef<-fit@fit$coefficients[,!colnames(fit@fit$coefficients) %in% c("(Intercept)", "scalingFactor")]
    pval<-fit@eb$p.value[,!colnames(fit@fit$coefficients) %in% c("(Intercept)", "scalingFactor")]
    coef.vector<-rename.features(coef, 'coef')
    pvalue.vector<-rename.features(pval, 'pval')
    paras<-cbind.data.frame(coef.vector, pvalue.vector)
    paras<-paras[, !duplicated(colnames(paras))]
  } 
  else{
    coef<-fit@fit$coefficients[,!colnames(fit@fit$coefficients) %in% c("(Intercept)", "scalingFactor")]
    pval<-fit@eb$p.value[,!colnames(fit@fit$coefficients) %in% c("(Intercept)", "scalingFactor")]
    paras<-cbind.data.frame(coef,pval)
    colnames(paras)<-c('coef', 'pval')
    paras$feature<-rownames(paras)
    paras$metadata<- names(metadata)
  }
  
  ###############################################
  # Calculate multiple qvalues only if prompted #
  ###############################################
  
  if(multiple_qvalues){
    paras<-append_qvalues(features, metadata, paras)
  } else{
    paras$qval_BH<-as.numeric(p.adjust(paras$pval, method = 'BH'))
  }
  
  #################
  # Return output #
  #################
  
  paras<-paras[order(paras$qval_BH, decreasing=FALSE),]
  paras<-dplyr::select(paras, c('feature', 'metadata'), everything())
  rownames(paras)<-NULL
  return(paras)  
}


###########################################
# Fit metagenomeSeq To A List of Datasets #
###########################################

list.metagenomeSeq<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.metagenomeSeq", "rename.features"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        "reshape2", "metagenomeSeq"),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.metagenomeSeq(features, metadata, libSize, ID, transformation, multiple_qvalues)
      DD$pairwiseAssociation<-paste('pairwiseAssociation', 1:nrow(DD), sep='')
      wh.TP<-intersect(grep("[[:print:]]+\\_TP$", DD$metadata), grep("[[:print:]]+\\_TP$", DD$feature))
      newname<-paste0(DD$pairwiseAssociation[wh.TP], "_TP")
      DD$pairwiseAssociation[wh.TP]<-newname
      DD<-dplyr::select(DD, c('pairwiseAssociation', 'feature', 'metadata'), everything())
      stop.time<-Sys.time()
      time<-as.numeric(round(difftime(stop.time, start.time, units="min"), 3), units = "mins")
      DD$time<-time
      return(DD)
    }
}



