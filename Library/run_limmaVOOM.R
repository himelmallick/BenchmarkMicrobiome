#############
# limmaVOOM #
#############

###########################
# Load Essential Packages #
###########################

# pacman, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
load_essential_packages()

###########################################
# Load Dedicated Method-specific Packages #
###########################################

pacman::p_load('reshape2')

if(! require("limma")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("limma")
}
suppressPackageStartupMessages(library(limma))

##############################
# Fit limmaVOOM To A Dataset #
##############################

fit.limmaVOOM = function(features, 
                         metadata, 
                         libSize, 
                         ID, 
                         transformation,
                         multiple_qvalues) {
  
  #########################
  # Transformation if any #
  #########################
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default limmaVOOM model. Use NONE.')
  
  ###############################
  # Standard limmaVOOM pipeline #
  ###############################
  
  x<-t(as.matrix(features)+1) # Convert to matrix, round up to nearest integer, and transpose
  design <- model.matrix(~., data=metadata)
  y <- voom(x,design,plot=FALSE)
  
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

#######################################
# Fit limmaVOOM To A List of Datasets #
#######################################

list.limmaVOOM<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.limmaVOOM", "rename.features"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        "reshape2", "limma"),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.limmaVOOM(features, metadata, libSize, ID, transformation, multiple_qvalues)
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

