##########
# DESeq2 #
##########

###########################
# Load Essential Packages #
###########################

# pacman, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
load_essential_packages()

###########################################
# Load Dedicated Method-specific Packages #
###########################################

pacman::p_load('reshape2')

if(! require("DESeq2")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
suppressPackageStartupMessages(library(DESeq2))

###########################
# Fit DESeq2 To A Dataset #
###########################

fit.DESeq2<-function(features, 
                     metadata, 
                     libSize, 
                     ID, 
                     transformation,
                     multiple_qvalues){
  
  #########################
  # Transformation if any #
  #########################
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default DESeq2 model. Use NONE.')

  ############################
  # Standard DESeq2 pipeline #
  ############################
  
  formula <- as.formula(paste('~', paste(colnames(metadata), collapse = "+"), sep=''))
  x <- DESeqDataSetFromMatrix(countData = t(as.matrix(features)), colData = metadata, design = formula)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(x), 1, gm_mean)
  x = estimateSizeFactors(x, geoMeans = geoMeans)
  fit <- DESeq(x)
  
  ###################
  # Combine results #
  ###################
  
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(coef(fit)[,-1], 'coef')
    pvalMatrix<-get_pval_DESeq2(fit)
    pvalue.vector<-rename.features(pvalMatrix[,-1], 'pval')
    paras<-cbind.data.frame(coef.vector, pvalue.vector)
    paras<-paras[, !duplicated(colnames(paras))]
    }
  else{
    coef<-coef(fit)[,-1]
    pval<-results(fit,name=resultsNames(fit)[2])$pvalue
    paras<-cbind.data.frame(coef, pval)
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
  
  paras<-paras[order(paras$qval_BH, decreasing = FALSE),]
  paras<-dplyr::select(paras, c('feature', 'metadata'), everything())
  rownames(paras)<-NULL
  return(paras)   
}

####################################
# Fit DESeq2 To A List of Datasets #
####################################

list.DESeq2<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.DESeq2", "rename.features", "get_pval_DESeq2"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        "reshape2", "DESeq2"),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.DESeq2(features, metadata, libSize, ID, transformation, multiple_qvalues)
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

