##########
# ALDEx2 #
##########

###########################
# Load Essential Packages #
###########################

# pacman, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
load_essential_packages()

###########################################
# Load Dedicated Method-specific Packages #
###########################################

if(! require("ALDEx2")) {
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install("ALDEx2")
  }

##########################
# Fit ANCOM To A Dataset #
##########################

fit.ALDEx2 = function(features, 
                      metadata, 
                      libSize, 
                      ID, 
                      transformation,
                      multiple_qvalues) {
  
  #########################
  # Transformation if any #
  #########################

  if (transformation!='NONE') stop ('Transformation currently not supported for a default ANCOM model. Use NONE.')

  ########################################################
  # ALDEx2 standard pipeline from Hawinkel et al. (2019) #
  ########################################################
  
  name_metadata <- names(metadata)
  grp <- metadata[,name_metadata]
  data<- as.data.frame(t(features))
  res = aldex(data, conditions = grp, test ="t")
  
  ###################
  # Combine results #
  ###################
  
  paras<-cbind.data.frame(coef = res$effect, pval = res$we.ep)
  paras$feature<-rownames(res)
  paras$metadata<- names(metadata)
  
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

####################################
# Fit ALDEx2 To A List of Datasets #
####################################

list.ALDEx2<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.ALDEx2"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        'ALDEx2'),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.ALDEx2(features, metadata, libSize, ID, transformation, multiple_qvalues)
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


