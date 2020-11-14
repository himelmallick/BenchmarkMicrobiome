##########
# ANCOM2 #
##########

###########################
# Load Essential Packages #
###########################

# pacman, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
load_essential_packages()

###########################################
# Load Dedicated Method-specific Packages #
###########################################

pacman::p_load('exactRankTests', 'openxlsx', 'DT', 'coin', 'compositions') # tidyverse is deliberastely excluded
devtools::source_url("https://github.com/FrederickHuangLin/ANCOM/blob/master/scripts/ancom_v2.1.R?raw=TRUE")

##########################
# Fit ANCOM To A Dataset #
##########################

fit.ANCOM2 = function(features, 
                      metadata, 
                      libSize, 
                      ID, 
                      transformation,
                      multiple_qvalues) {
  
  #########################
  # Transformation if any #
  #########################

  if (transformation!='NONE') stop ('Transformation currently not supported for a default ANCOM model. Use NONE.')

  #############################################################################################
  # ANCOM standard pipeline for DA (adopted from https://github.com/FrederickHuangLin/ANCOM) #
  #############################################################################################
  
  ###############################
  # Step 1: ANCOM preprocessing #
  ###############################

  features_ancom<-as.data.frame(t(features))
  metadata_ancom<-rownames_to_column(metadata, 'ID')
  preprocess.ancom = feature_table_pre_process(feature_table = features_ancom, 
                                     meta_data = metadata_ancom, 
                                     sample_var = "ID", 
                                     group_var = NULL,
                                     out_cut = 0.05,
                                     zero_cut = 0.90,
                                     lib_cut = 1000,
                                     neg_lb = FALSE)
  feature_table = preprocess.ancom$feature_table 
  meta_data = preprocess.ancom$meta_data 
  struc_zero = preprocess.ancom$structure_zeros 

  
  #####################
  # Step 2: Run ANCOM #
  #####################

  if(!length(ID)==length(unique(ID))){
    res <- ANCOM(feature_table = feature_table, 
                 meta_data = meta_data, 
                 struc_zero = struc_zero,
                 main = colnames(metadata),
                 alpha = 0.05,
                 adj_formula = NULL,
                 rand_formula = "~ 1 | ID",
                 p_adj_method = 'BH') 
  } else{
    res <- ANCOM(feature_table = feature_table, 
                 meta_data = meta_data, 
                 struc_zero = struc_zero,
                 main = colnames(metadata),
                 alpha = 0.05,
                 adj_formula = NULL,
                 rand_formula = NULL,
                 p_adj_method = 'BH')
    }

  ###################
  # Combine Results #
  ###################
  
  paras <- data.frame(coef = res$out)
  paras$pval<-1 # Fake p-values
  paras$feature = colnames(features)
  paras$metadata = names(metadata)
  paras$pval[paras$coef.detected_0.7] <- 0 # Fake p-values
  
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
# Fit ANCOM2 To A List of Datasets #
####################################

list.ANCOM2<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE)){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.ANCOM2"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        'exactRankTests', 'openxlsx', 'DT', 'coin', 'compositions'),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.ANCOM2(features, metadata, libSize, ID, transformation, multiple_qvalues)
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


