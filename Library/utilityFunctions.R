###########################
# Load Essential Packages #
###########################

load_essential_packages<-function(){
  
  if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
  suppressPackageStartupMessages(library("pacman"))
  pacman::p_load("devtools")
  
  ###################
  # Core R Packages #
  ###################
  
  pacman::p_load('tidyverse', 'fdrtool', 'ashr')
  
  ########################
  # Core GitHub Packages #
  ########################
  
  if(! require("GMPR")) {
    devtools::install_github("lichen-lab/GMPR")
  }
  if(! require("swfdr")) {
    devtools::install_github("leekgroup/swfdr")
  }
  library(GMPR)
  library(swfdr)
  
  ##############################
  # Core Bioconductor Packages #
  ##############################
  
  if(! require("genefilter")) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("genefilter")
  }
  library(genefilter)
  
  if(! require("IHW")) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("IHW")
  }
  library(IHW)
}

####################################################################################
# Given features, metadata, and summary table consisting of coef, stderr, p-values #
# from a given model append a range of qvalues as new columns to the summary table #
####################################################################################

append_qvalues<-function(features, metadata, paras){
  
  ########################
  # Define all functions #
  ########################
  
  try_my_func <- function() {
    tryCatch(my_func(), error = function(err){NA})
  }
  
  ##############################################################
  # Calculate normalized means and independent filtering index #
  ##############################################################
  
  df_ashr<-nrow(metadata) - ncol(metadata) - 1 # DF for ashr (t)
  gmpr.size.factor<-GMPR(features)
  features_norm <- features/gmpr.size.factor
  baseMean<-colMeans(features_norm)
  qval_IF<-as.numeric(pvalueAdjustment_HM(baseMean = baseMean, 
                                          pValue = paras$pval, 
                                          pAdjustMethod = 'BH')$padj) # Independent filtering index
  
  ##########################################################
  # Calculate multiple FDR-adjusted p-values when possible #
  ##########################################################
  
    paras<- paras %>% 
    dplyr::mutate(qval_BH = tryCatch(p.adjust(pval, method = 'BH'),
                                     error = function(err){NA})) %>% # Benjamini-Hochberg (BH)
    
    dplyr::mutate(qval_BY = tryCatch(p.adjust(pval, method = 'BY'),
                                     error = function(err){NA})) %>% # Benjamini-Yekutieli (BY)
    
    dplyr::mutate(qval_LFDR = tryCatch(fdrtool::fdrtool(pval, 
                                               statistic="pvalue")$lfdr,
                                       error = function(err){NA})) %>%  # Storey's local FDR (LFDR)
    
    dplyr::mutate(qval_ST = tryCatch(qvalue::qvalue(pval)$qvalues,
                                     error = function(err){NA})) %>%  # Storey's q-value (ST)
    
    dplyr::mutate(qval_BL = tryCatch(swfdr::lm_qvalue(pval, 
                                             X = baseMean)$qvalues),
                                      error = function(err){NA}) %>% # Boca and Leek (BL)
    
    dplyr::mutate(qval_IHW = tryCatch(IHW::ihw(pvalues = pval, 
                                      covariates = IHW::groups_by_filter(baseMean, 20), 
                                      alpha =  1)@df$adj_pvalue),
                                      error = function(err){NA})  %>% # IHW q-value (IHW)
    
    dplyr::mutate(qval_ASH1 = tryCatch(ashr::get_qvalue(ashr::ash(coef, stderr)),
                                       error = function(err){NA})) %>% # ashr q-value (ASH1)
    
    dplyr::mutate(qval_ASH2 = tryCatch(ashr::get_svalue(ashr::ash(coef, stderr)),
                                       error = function(err){NA})) %>% # ashr s-value (ASH2)
    
    dplyr::mutate(qval_ASH3 = tryCatch(ashr::get_svalue(ashr::ash(coef, 
                                                                  stderr, 
                                                                  method = 'shrink')),
                                       error = function(err){NA})) %>% # ashr regularized s-value (ASH3)
    
    dplyr::mutate(qval_ASH4 = tryCatch(ashr::get_svalue(ashr::ash(coef, 
                                                                  stderr, 
                                                                  method = 'shrink', 
                                                                  df = df_ashr)),
                                       error = function(err){NA})) %>% # ashr t-regularized s-value (ASH4)
  
    
  ##############################################################################
  # Replace with NA if NA if detected as outlier by Independent Filtering (IF) #
  ##############################################################################
  
    dplyr::mutate(qval_BH.IF =  tryCatch(ifelse(is.na(qval_IF), NA, qval_BH),
                                         error = function(err){NA})) %>% # BH
    
    dplyr::mutate(qval_BY.IF =  tryCatch(ifelse(is.na(qval_IF), NA, qval_BY),
                                         error = function(err){NA})) %>% # BY
    
    dplyr::mutate(qval_LFDR.IF =  tryCatch(ifelse(is.na(qval_IF), NA, qval_LFDR),
                                           error = function(err){NA})) %>% # LFDR
    
    dplyr::mutate(qval_ST.IF =  tryCatch(ifelse(is.na(qval_IF), NA, qval_ST),
                                         error = function(err){NA})) %>% # ST
    
    dplyr::mutate(qval_BL.IF =  tryCatch(ifelse(is.na(qval_IF), NA, qval_BL),
                                         error = function(err){NA})) %>% # BL
    
    dplyr::mutate(qval_IHW.IF =  tryCatch(ifelse(is.na(qval_IF), NA, qval_IHW),
                                          error = function(err){NA})) %>% # IHW
    
    dplyr::mutate(qval_ASH1.IF =  tryCatch(ifelse(is.na(qval_IF), NA, qval_ASH1),
                                           error = function(err){NA})) %>% # ASH1
    
    dplyr::mutate(qval_ASH2.IF =  tryCatch(ifelse(is.na(qval_IF), NA, qval_ASH2),
                                           error = function(err){NA})) %>% # ASH2
    
    dplyr::mutate(qval_ASH3.IF =  tryCatch(ifelse(is.na(qval_IF), NA, qval_ASH3),
                                           error = function(err){NA})) %>% # ASH3
    
    dplyr::mutate(qval_ASH4.IF =  tryCatch(ifelse(is.na(qval_IF), NA, qval_ASH4),
                                           error = function(err){NA}))     # ASH4

  return(paras)

}

#############################################
# Helper function for Independent Filtering #
#############################################

pvalueAdjustment_HM <- function(baseMean, 
                                filter, 
                                pValue,
                                theta, 
                                alpha = 0.05, 
                                pAdjustMethod = "BH") {
  
  ##############################################################
  # Adapted from the pValueAdjustment function from the DESeq2 #
  ##############################################################
  
  #######################
  # Extract user inputs #
  #######################
  
  if (missing(filter)) {
    filter <- baseMean
  }
  if (missing(theta)) {
    lowerQuantile <- mean(filter == 0)
    if (lowerQuantile < .95) upperQuantile <- .95 else upperQuantile <- 1
    theta <- seq(lowerQuantile, upperQuantile, length=50)
  }
  
  #################################
  # Perform independent filtering #
  #################################
  
  stopifnot(length(theta) > 1)
  filtPadj <- filtered_p(filter=filter, test=pValue,
                         theta=theta, method=pAdjustMethod)
  numRej  <- colSums(filtPadj < alpha, na.rm = TRUE)
  # prevent over-aggressive filtering when all genes are null,
  # by requiring the max number of rejections is above a fitted curve.
  # If the max number of rejection is not greater than 10, then don't
  # perform independent filtering at all.
  lo.fit <- lowess(numRej ~ theta, f=1/5)
  if (max(numRej) <= 10) {
    j <- 1
  } else {
    residual <- if (all(numRej==0)) {
      0
    } else {
      numRej[numRej > 0] - lo.fit$y[numRej > 0]
    }
    thresh <- max(lo.fit$y) - sqrt(mean(residual^2))
    j <- if (any(numRej > thresh)) {
      which(numRej > thresh)[1]
    } else {
      1
    }
  }
  padj <- filtPadj[, j, drop=TRUE]
  cutoffs <- quantile(filter, theta)
  filterThreshold <- cutoffs[j]
  filterNumRej <- data.frame(theta=theta, numRej=numRej)
  filterTheta <- theta[j]
  
  return(list(padj=padj, filterThreshold=filterThreshold, filterTheta=filterTheta, filterNumRej = filterNumRej, lo.fit=lo.fit, alpha=alpha))
  
}

######################
## RFA Normalization #
######################

list.RFA<- function(physeq){
  foreach(physeq = physeq, .packages = "phyloseq") %dopar% {
    physeq$features = rarefy_even_depth(otu_table(physeq$features, taxa_are_rows = FALSE))
  return(physeq)
  }
}

######################
## CSS Normalization #
######################

# Apply CSS Normalization To A Dataset

CSSnorm = function(physeq) {
  
  # Extract ID and Metadata and Keep Them Unchanged
  ID<-physeq$ID
  metadata<-physeq$metadata
  
  # Extract Features
  features = as.matrix(physeq$features)
  dd<-colnames(features)
  
  # CSS Normalizing the Data
  # Create the metagenomeSeq object
  MGS = newMRexperiment(t(features), featureData=NULL, libSize=NULL, normFactors=NULL)
  # Trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
  MGS = cumNorm(MGS, p=cumNormStat(MGS))
  # Calculate scaling factors
  libSize <- normFactors(MGS)
  # Save the normalized data as data.frame
  features_norm = as.data.frame(t(MRcounts(MGS, norm=TRUE, log=FALSE))) 
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_norm) <- dd;
  
  # Return as list
  return(list(metadata=metadata, features=features_norm, ID=ID, libSize=libSize))
}

# Apply CSS Normalization To A List of Datasets

list.CSS<-function(physeq){
  foreach(physeq = physeq, .packages = "metagenomeSeq", .export="CSSnorm", .errorhandling = 'remove') %dopar% {
    return(CSSnorm(physeq))
  }
}

######################
## TSS Normalization #
######################

# Apply TSS Normalization To A Dataset

TSSnorm = function(physeq) {
  
  # Extract ID and Metadata and Keep Them Unchanged
  ID<-physeq$ID
  metadata<-physeq$metadata
  
  # Extract Features
  features = as.matrix(physeq$features)
  dd<-colnames(features)
  
  # TSS Normalizing the Data
  features_norm <- decostand(features, method="total", MARGIN=1)
  
  # Convert back to data frame
  features_norm<-as.data.frame(features_norm)
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_norm) <- dd;
  
  # Force the Sequencing Depth to be the Original Library Size (Same as CLR)
  libSize<-rowSums(features)
  
  # Return as list
  return(list(metadata=metadata, features=features_norm, ID=ID, libSize=libSize))
}

# Apply TSS Normalization To A List of Datasets

list.TSS<-function(physeq){
  foreach(physeq = physeq, .packages = "vegan", .export="TSSnorm", .errorhandling = 'remove') %dopar% {
    return(TSSnorm(physeq))
  }
}

######################
## CLR Normalization #
######################

# Apply CLR Normalization To A Dataset

CLRnorm = function(physeq) {
  
  # Extract ID and Metadata and Keep Them Unchanged
  ID<-physeq$ID
  metadata<-physeq$metadata
  
  # Extract Features
  features = as.matrix(physeq$features)
  dd<-colnames(features)
  
  # Force the Sequencing Depth to be the Original Library Size (Same as TSS)
  libSize<-rowSums(features)

  # CLR
  features_norm_clr<-clr(features+1)
  
  # Convert back to data frame
  features_norm_clr<-as.data.frame(features_norm_clr)
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_norm_clr) <- dd;
  
  # Return as list
  return(list(metadata=metadata, features= features_norm_clr, ID=ID, libSize=libSize))
}

# Apply CLR Normalization To A List of Datasets

list.CLR<-function(physeq){
  foreach(physeq = physeq, 
          .packages = c("vegan", "chemometrics"),
          .export="CLRnorm", .errorhandling = 'remove') %dopar% {
            return(CLRnorm(physeq))
          }
}

######################
## TMM Normalization #
######################

# Apply TMM Normalization To A Dataset

TMMnorm = function(physeq) {
  
  # Extract ID and Metadata and Keep Them Unchanged
  ID<-physeq$ID
  metadata<-physeq$metadata
  
  # Extract Features
  features = as.matrix(physeq$features)
  dd<-colnames(features)
  
  # TMM Normalizing the Data
  X<-t(features);
  libSize = edgeR::calcNormFactors(X,method="TMM") #Calculate normaization factors
  eff.lib.size = colSums(X)*libSize;
  ref.lib.size = mean(eff.lib.size); #Use the mean of the effective library sizes as a reference library size 
  X.output = sweep(X,MARGIN=2,eff.lib.size,"/")*ref.lib.size; #Normalized read counts
  
  # Convert back to data frame
  features_norm<-as.data.frame(t(X.output))
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_norm) <- dd;
  
  # Return as list
  return(list(metadata=metadata, features=features_norm, ID=ID, libSize=libSize))
}

# Apply TMMnorm Normalization To A List of Datasets

list.TMM<-function(physeq){
  foreach(physeq = physeq, .packages = "edgeR", .export="TMMnorm", .errorhandling = 'remove') %dopar% {
    return(TMMnorm(physeq))
  }
}

######################
## RLE Normalization #
######################

# Apply RLE Normalization To A Dataset

RLEnorm = function(physeq) {
  
  # Extract ID and Metadata and Keep Them Unchanged
  ID<-physeq$ID
  metadata<-physeq$metadata
  
  # Extract Features
  features = as.matrix(physeq$features)
  dd<-colnames(features)
  
  # RLE Normalizing the Data
  formula <- as.formula(paste('~', paste(colnames(metadata), collapse = "+"), sep=''))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = t(features), colData = metadata, design = formula)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(dds), 1, gm_mean)
  dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  X.output<-counts(dds, normalized=TRUE)
  
  # Convert back to data frame
  features_norm<-as.data.frame(t(X.output))
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_norm) <- dd;
  
  # Extract Scaling Factors
  libSize<-sizeFactors(dds)
  
  # Return as list
  return(list(metadata=metadata, features=features_norm, ID=ID, libSize=libSize))
}

# Apply RLE Normalization To A List of Datasets

list.RLE<-function(physeq){
  foreach(physeq = physeq, .packages = "DESeq2", .export="RLEnorm", .errorhandling = 'remove') %dopar% {
    return(RLEnorm(physeq))
  }
}


######################
# GMPR Normalization #
######################

# Apply GMPR Normalization To A Dataset

GMPRnorm = function(physeq) {
  
  # Extract ID and Metadata and Keep Them Unchanged
  ID<-physeq$ID
  metadata<-physeq$metadata
  
  # Extract Features
  features = as.matrix(physeq$features)
  
  # GMPR Normalizing the Data
  gmpr.size.factor<-GMPR::GMPR(features)
  
  # Reset library size
  libSize<-gmpr.size.factor
  
  # Return as list
  return(list(metadata=metadata, features=features, ID=ID, libSize=libSize))
}

# Apply GMPR Normalization To A List of Datasets

list.GMPR <-function(physeq){
  foreach(physeq = physeq, .packages = "GMPR", .export="GMPRnorm", .errorhandling = 'remove') %dopar% {
    return(GMPRnorm(physeq))
  }
}

#######################
# SCRAN Normalization #
#######################

# Apply SCRAN Normalization To A Dataset

SCRANnorm = function(physeq) {
  
  # Extract ID and Metadata and Keep Them Unchanged
  ID<-physeq$ID
  metadata<-physeq$metadata
  
  # Extract Features
  features = as.matrix(physeq$features)
  
  # SCRAN Normalizing the Data
  sce<-t(features)
  clusters <- scran::quickCluster(sce)
  sizeFactors <- scran::computeSumFactors(sce, clusters=clusters)
  
  # Reset library size
  libSize<-sizeFactors
  
  # Return as list
  return(list(metadata=metadata, features=features, ID=ID, libSize=libSize))
}

# Apply SCRAN Normalization To A List of Datasets

list.SCRAN <-function(physeq){
  foreach(physeq = physeq, .packages = "scran", .export="SCRANnorm", .errorhandling = 'remove') %dopar% {
    return(SCRANnorm(physeq))
  }
}


#########################
# Summarization Scripts #
#########################

################################################################################
# Calculate evaluation metrics for a single dataset but multiple (alpha, qval) #
################################################################################

eval_res_list_repeated<-function(resi, 
                                 qval_Type = c('qval_BH', 'qval_BY', 'qval_LFDR',
                                               'qval_ST', 'qval_BL', 'qval_IHW',
                                               'qval_ASH1', 'qval_ASH2', 
                                               'qval_ASH3', 'qval_ASH4',
                                               'qval_BH.IF', 'qval_BY.IF', 
                                               'qval_LFDR.IF','qval_ST.IF', 
                                               'qval_BL.IF', 'qval_IHW.IF',
                                               'qval_ASH1.IF', 'qval_ASH2.IF', 
                                               'qval_ASH3.IF', 'qval_ASH4.IF'),
                                 alpha_range = c(0, 0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.50, 0.75, 1.00)) {
  
  
  ########################
  # Intialize data frame #
  ########################
  
  df_final = data.frame()
  
  ###################################################################
  # Calculate evaluation metrics for each (alpha, qval) combination #
  ###################################################################
  
  for (i in seq_along(qval_Type)){
    for (j in seq_along(alpha_range)){ 
      temp<-eval_res_list(resi, qval = qval_Type[i], alpha = alpha_range[j])
      df <- data.frame(t(temp))
      df_final <- data.frame(rbind(df_final,df))
    }
    }
  
  ##########
  # Return #
  ##########
  
  return(df_final)
}

#######################################################################
# Calculate evaluation metrics for a single dataset and (alpha, qval) #
#######################################################################

eval_res_list = function(resi, qval = "qval_BH", alpha = 0.05) {
  
  # Define a q-value column
  names(resi)[names(resi) == qval] <- 'qval'
  
  # Replace Missing Q-values to the Highest Possible Value 1.0
  resi[is.na(resi[, "qval"]), "qval"] <- 1
  
  # Evaluate Detection Performance
  time = mean(resi[,"time"], na.rm=TRUE)
  wh.pred = (resi[, "qval"] <= alpha)
  wh.pos = which(wh.pred)
  wh.neg = which(!wh.pred)
  wh.TP = grep("[[:print:]]+\\_TP$", resi[, "pairwiseAssociation"])
  FPs = sum(!wh.pos %in% wh.TP)
  TPs = sum(wh.pos %in% wh.TP)
  TNs = sum(!wh.neg %in% wh.TP)
  FNs = sum(wh.neg %in% wh.TP)
  
  
  # Sensitivity: True Positives Divided by All Positives (Sum of True
  # Positives and False Negatives)
  Sensitivity = TPs/(TPs + FNs)
  
  # Specificity: True Negatives Divided by All Negatives (Sum of True
  # Negatives and False Positives)
  Specificity = TNs/(TNs + FPs)
  
  # False Discovery Rate: False Positives Divided by All Detected Positives
  FDR = if ((TPs + FPs) == 0) 
    0 else FPs/(TPs + FPs)
  # If no true positives, return NA's for irrelevant measures
  
  # FScore
  FScore<-2*TPs/(2*TPs + FNs + FPs) 
  
  # Matthew's Correlation Coefficient  
  numerator <- (TPs * TNs - FPs * FNs)
  denominator <- sqrt(as.numeric((TPs + FPs))*as.numeric((TPs + FNs))*as.numeric((TNs + FPs))*as.numeric((TNs + FNs)))
  if(denominator == 0) denominator <- 1
  if(is.na(denominator)) denominator <- 1
  MCC <- numerator/denominator
  
  # AUC, pAUC (FPR < 0.20)
  wh.truth = (1:nrow(resi) %in% wh.TP)
  if (all(wh.truth=='TRUE')) {
    AUC=1
    pAUC=1
    fAUC=1
  } else if (all(wh.truth=='FALSE')) {
    AUC=0
    pAUC=0
  } else {
    pred <- prediction(as.numeric(wh.pred), factor(wh.truth, levels=c("TRUE", "FALSE")))
    AUC = performance(pred, "auc")@y.values[[1]]
    pAUC = performance(pred, "auc", fpr.stop=0.20)@y.values[[1]][[1]]
  }
  
  # Departure from Uniformity Under the Null
  AucAocVals <- AucAocFun(resi[!wh.truth, "pval"], plotIt = FALSE, pch = 20, type = "l")
  totalArea = AucAocVals["conservArea"] + AucAocVals["liberalArea"]
  names(totalArea) = "totalArea"  
  
  # Return
  return(c(Sensitivity = Sensitivity,
           Specificity = Specificity,
           FDR = FDR,
           FScore = FScore,
           MCC = MCC,
           AUC = AUC,
           pAUC = pAUC,
           AucAocVals["conservArea"],
           AucAocVals["liberalArea"],
           totalArea,
           time = time,
           alpha = alpha,
           qval = sub(".*_", "", qval)))
}

#############################################
# Extract Parameters from A List of Results #
#############################################

make_power_df = function(reslist, simparamslabels) {
  powerdf = ldply(reslist, .parallel=TRUE)
  colnames(powerdf)[1] <- "Combinations"
  paramdf = ldply(strsplit(powerdf[, "Combinations"], "_"), .parallel=TRUE)
  colnames(paramdf) <- simparamslabels
  powerdf = cbind(powerdf, paramdf)
  return(powerdf)
}

#################################################
# Organize All Lists into A Coherent Data Frame #
#################################################

list.perfdf<-function(resultslist, simparams, simparamslabels){
  foreach(resultslist = resultslist, .export=c("make_power_df", "eval_res_list", "AucAocFun"),
          .packages = c("ROCR", "plyr"), .errorhandling = 'remove') %dopar% {
  perflist <- lapply(resultslist, eval_res_list)
  if (is.null(names(resultslist))) {
    names(perflist) <- simparams
  } else {
    names(perflist) <- names(resultslist)
  }
  perfdf = make_power_df(perflist, simparamslabels)
  return(perfdf)
}}


#######################################
# Arc-Sine Square Root Transformation #
#######################################

# Arc Sine Square Root Transformation
AST<-function(x){
  return(sign(x)*asin(sqrt(abs(x))))
}

########################
# Logit Transformation #
########################

# LOGIT Transformation 
LOGIT<-function(x){
  y<-car::logit(x, adjust=0)
  y[!is.finite(y)]<-0
  return(y)
}

# Shifted Logit Transformation (Lukens et al, 2014, Nature)
# LOGIT_S<-function(x){
#   y<-0.5*log(x/(1-x)) + 10
#   y[!is.finite(y)]<-0
#   return(y)
# }

######################
# Log Transformation #
######################

# LOG Transformation
LOG<-function(x){
  return(log(x+1))
}

########################################################################
# Apply SimpleFilter (Prevalence and Abundance Filtering) To A Dataset #
########################################################################

SimpleFilter = function(physeq, Threshold_Abundance, Threshold_Prevalence) {
  features<-physeq$features
  filtered_features<-features[,colSums(features > Threshold_Abundance) > nrow(features)*Threshold_Prevalence] 
  physeq$features<-filtered_features
  return(physeq)
}

#################################################################################
# Apply SimpleFilter (Prevalence and Abundance Filtering) To A List of Datasets #
#################################################################################

list.SimpleFilter<- function(simlist, Threshold_Abundance, Threshold_Prevalence){
  return(lapply(simlist, SimpleFilter, Threshold_Abundance, Threshold_Prevalence))
}

#############################################
# Calculate No. of TP Features in A Dataset #
#############################################

TruePositives<-function(physeq, breakitdown) {
  features<-physeq$features
  d<-colnames(features)[grep("[[:print:]]+\\_TP$", colnames(features))]
  if(breakitdown==TRUE) {
    return(d)
  } else {
    return(length(d))
  }
}

######################################################
# Calculate No. of TP Features in A List of Datasets #
#####################################################

list.TruePositives<- function(simlist, breakitdown){
  return(lapply(simlist, TruePositives, breakitdown))
}

#####################################################
# Calculate Prevalence of All Features in A Dataset #
#####################################################

calculatePrevalence<-function(physeq) {
  features<-physeq$features
  d<-colSums(features > Threshold_Abundance)/nrow(features)
  return(d)
}

##############################################################
# Calculate Prevalence of All Features in A List of Datasets #
##############################################################

list.calculatePrevalence<- function(simlist){
  return(lapply(simlist, calculatePrevalence))
}

#################################
# Check if a metadata is binary #
#################################

is.binary <- function(v) {
  x <- unique(v)
  dim(x)[1] - sum(is.na(x)) == 2L
}

########################################
# Function from Hawinkel et al. (2017) #
########################################

######################################################################################
# A function to look at the distribution of the P-values under the null distribution #
######################################################################################

### Compute AUC and AOC of P-value distribution to be applied on the Monte
### Carlo distribution of each feature separately
AucAocFun <- function(pVals, maxVal = 0.25, plotIt = FALSE, ...) {
  ## sum of height differences between the curve and the y=x line, half maximum
  ## = 0.5 * max(pVals) * length(pVals) this is the case of 100% rejection with
  ## any alpha value _onlyTotArea_ controls if only the total area is returned,
  ## and not also the Area Over and Under the y=x line.  halfMaxArea <- 0.5 *
  ## max(maxVal) * length(pVals)
  halfMaxArea <- 0.5 * length(pVals)
  
  pVals <- pVals[!is.na(pVals)]
  estimated <- sort(pVals)
  theoretic <- seq_along(pVals)/length(pVals)
  
  if (plotIt) {
    plot(theoretic, estimated, ...)
    abline(0, 1)
  } else {
  }
  
  diffPVals <- theoretic - estimated
  indConserv <- theoretic <= estimated
  conservArea <- sum(-diffPVals[indConserv])/halfMaxArea
  liberalArea <- sum(diffPVals[!indConserv])/halfMaxArea
  
  c(conservArea = conservArea, liberalArea = liberalArea, totalArea = liberalArea + 
      conservArea)
  
}  # END - function: aucAocFun, AUC/AOC calculation for p-values distribution

###########################################
# Helper Functions for Individual Methods #
###########################################

########################################################
# Convert A Limma/edgeR Output Matrix into Long Format #
########################################################

rename.features<-function(x, name){
  x<-reshape2::melt(x);
  colnames(x)<-c('feature', 'metadata', name);
  return(x);
}

################################
# Get P-values from A MAST Fit #
################################

get_pval_MAST<-function(fit){
  mat<-as.matrix(fit[,3,1])
  rownames(mat)<- row.names(fit[,,1])
  colnames(mat)<- names(metadata)
  
  rownames<-rownames(mat)
  colnames<-colnames(mat)
  n<-dim(mat)[2]
  List <- list()
  for(i in 1:n){
    List[[i]] <- fit[,3,2+i]
  }
  Matrix = do.call(cbind, List)
  rownames(Matrix)<-rownames
  colnames(Matrix)<-colnames
  return(Matrix)
}

##################################
# Get P-values from An edgeR Fit #
##################################

get_pval_edgeR<-function(fit){
  mat<-as.matrix(fit$coefficients)
  rownames<-rownames(mat)
  colnames<-colnames(mat)
  n<-dim(fit$coefficients)[2]
  List <- list()
  for(i in 1:n){
    List[[i]] <- glmLRT(fit, i)$table$PValue
  }
  Matrix = do.call(cbind, List)
  rownames(Matrix)<-rownames
  colnames(Matrix)<-colnames
  return(Matrix)
}

##################################
# Get P-values from A DESeq2 Fit #
##################################

get_pval_DESeq2<-function(fit){
  List <- list()
  for(i in 1:length(resultsNames(fit))){
    List[[i]] <- results(fit,name=resultsNames(fit)[i])$pvalue
  }
  Matrix = do.call(cbind, List)
  rownames(Matrix)<-names(fit)
  colnames(Matrix)<-resultsNames(fit)
  return(Matrix)
}
