############
# Spearman #
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
# Fit Spearman To A Dataset #
#############################

fit.Spearman  = function(features, 
                         metadata, 
                         libSize, 
                         ID, 
                         transformation) {
  
  #########################
  # Transformation if any #
  #########################
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default Spearman model. Use NONE.')

  ##############################
  # Standard Spearman pipeline #
  ##############################
  
  spe <- function(x){
    tryCatch({
    fit1 <- cor.test(x, as.vector(t(metadata)), method='spearman')
  }, error=function(err){
    fit1 <- try({cor.test(x, as.vector(t(metadata)), method='spearman')}) 
    return(fit1)
  })}
  
  #####################
  # Summarize Results #
  #####################
  
  spes <- apply(features, 2, spe)
  paras <- data.frame(pval = sapply(spes, function(x) x$p.value))
  paras$coef <- sapply(spes, function(x) x$statistic)
  paras$metadata<-colnames(metadata)
  paras$feature<-colnames(features)
  
  ###################
  # Combine results #
  ###################
  
  paras<-append_qvalues(features, metadata, paras)
  
  #################
  # Return output #
  #################
  
  paras<-paras[order(paras$qval_BH, decreasing=FALSE),]
  paras<-dplyr::select(paras, c('feature', 'metadata'), everything())
  rownames(paras)<-NULL
  return(paras)   
}


######################################
# Fit ANCOM2 To A List of Datasets #
######################################

list.ANCOM2<-function(physeq, transformation = 'NONE'){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.ANCOM2"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        'tibble', 'exactRankTests', 'openxlsx', 'DT', 'dplyr', 'coin', 'compositions'),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.ANCOM2(features, metadata, libSize, ID, transformation)
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