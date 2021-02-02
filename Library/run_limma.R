#########
# limma #
#########

###########################
# Load Essential Packages #
###########################

# pacman, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
load_essential_packages()

###########################################
# Load Dedicated Method-specific Packages #
###########################################

pacman::p_load('reshape2', 'car')

if(! require("limma")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("limma")
}
suppressPackageStartupMessages(library(limma))

##########################
# Fit limma To A Dataset #
##########################

fit.limma = function(features, 
                     metadata, 
                     libSize, 
                     ID, 
                     transformation,
                     multiple_qvalues) {
  
  #########################
  # Transformation if any #
  #########################
  
  if (!transformation %in% c('AST', 'LOG', 'LOGIT', 'NONE')) {
    stop ('Transformation should be one of AST/LOG/LOGIT/NONE for limma.')
  }
  
  if (transformation =='LOG')   {
    x<-apply(features, 2, LOG)
  }
  
  if (transformation =='LOGIT')   {
    x<-apply(features, 2, LOGIT)
  }
  
  if (transformation =='AST')   {
    x<-apply(features, 2, AST)
  }
  
  if (transformation =='NONE')   {
    x<-features
  }
  
  # Convert to matrix, and transpose
  y<-t(as.matrix(x)) 
  
  ###########################
  # Standard limma pipeline #
  ###########################
  
  design <- model.matrix(~., data = metadata)
  
  #############
  # Fit limma # 
  #############
  
  ############################
  # Random Effect Adjustment #
  ############################
  
  if (!length(ID)==length(unique(ID))){
    dupcor <-  limma::duplicateCorrelation(y, design, block = ID)
    fit <- limma::lmFit(y,design, block = ID, correlation = dupcor$cor)
  } else{ 
    fit <- limma::lmFit(y,design)}
  
  ##############################
  # Empirical Bayes Adjustment #
  ##############################
  
  fit <- limma::eBayes(fit)
  
  ###################
  # Combine results #
  ###################
  
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(fit$coefficients[,-1], 'coef')
    pvalue.vector<-rename.features(fit$p.value[,-1], 'pval')
    paras<-cbind.data.frame(coef.vector, pvalue.vector)
    paras<-paras[, !duplicated(colnames(paras))]
  } 
  else{
    coef<-fit$coefficients[,-1]
    pval<-fit$p.value[,-1]
    paras<-cbind.data.frame(coef,pval)
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

###################################
# Fit limma To A List of Datasets #
###################################

list.limma<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.limma", "rename.features", "AST", "LOG", "LOGIT"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        "reshape2", "limma"),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.limma(features, metadata, libSize, ID, transformation, multiple_qvalues)
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

