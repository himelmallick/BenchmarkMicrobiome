############
# Wilcoxon #
############

###########################
# Load Essential Packages #
###########################

# pacman, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
load_essential_packages()

###########################################
# Load Dedicated Method-specific Packages #
###########################################

#############################
# Fit Wilcoxon To A Dataset #
#############################

fit.Wilcoxon  = function(features, 
                         metadata, 
                         libSize, 
                         ID, 
                         transformation,
                         multiple_qvalues) {
  
  #########################
  # Transformation if any #
  #########################
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default Wilcoxon model. Use NONE.')
  
  ##############################
  # Standard Wilcoxon pipeline #
  ##############################
  
  group0indexes <- which(metadata == 0)
  group1indexes <- which(metadata == 1)
  
  spe <- function(x){

    tryCatch({
      fit1 <- wilcox.test(x[group0indexes], x[group1indexes])
    }, error=function(err){
      fit1 <- try({wilcox.test(x[group0indexes], x[group1indexes])}) 
      return(fit1)
    })}
  
  ###################
  # Combine results #
  ###################
  
  spes <- apply(features, 2, spe)
  paras <- data.frame(pval = sapply(spes, function(x) x$p.value))
  paras$coef <- sapply(spes, function(x) x$statistic)
  paras$metadata<-colnames(metadata)
  paras$feature<-colnames(features)
  
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

######################################
# Fit Wilcoxon To A List of Datasets #
######################################

list.Wilcoxon<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.Wilcoxon"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW"),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.Wilcoxon(features, metadata, libSize, ID, transformation, multiple_qvalues)
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