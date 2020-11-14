#########################################################
# Vanilla Linear Model (LM) with Library Size as Offset #
#########################################################

###########################
# Load Essential Packages #
###########################

# pacman, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
load_essential_packages()

###########################################
# Load Dedicated Method-specific Packages #
###########################################

pacman::p_load('pbapply', 'car', 'nlme')

#######################
# Fit LM To A Dataset #
#######################

fit.LM <- function(features, 
                   metadata, 
                   libSize, 
                   ID, 
                   transformation,
                   multiple_qvalues){
  
  ##################
  # Transformation #
  ##################
  
  if (!transformation %in% c('AST', 'LOG', 'LOGIT', 'NONE')) {
    stop ('Transformation should be one of AST/LOG/LOGIT/NONE for LM')
  }
  
  if (transformation =='LOG')   {
    features<-apply(features, 2, LOG)
  }

  if (transformation =='LOGIT')   {
    features<-apply(features, 2, LOGIT)
  }
  
  if (transformation =='AST')   {
    features<-apply(features, 2, AST)
  }
  
  if (transformation =='NONE')   {
    features<-features
  }
  
  #####################
  # Per-feature model #
  #####################
  
  paras <- pbapply::pbsapply(1:ncol(features), simplify=FALSE, function(x){
    
    ###############################
    # Extract features one by one #
    ###############################
    
    featuresVector <- features[, x]

    #################################
    # Create per-feature input data #
    #################################
    
    dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata, libSize, ID)
    formula<-as.formula(paste("expr ~ ", paste(colnames(metadata), collapse= "+")))
    
    ##############################################
    # Automatic library size adjustment for GLMs #
    ##############################################
    
    if(length(unique(libSize)) > 1){ # To prevent offsetting with TSS-normalized data 
      if (transformation =='LOG'){
        formula<-update(formula, . ~ . - offset(log(libSize)))
      } else{
        formula<-update(formula, . ~ . - offset(libSize))
      }
    }
  
    ########################
    # Random effects model #
    ########################
    
      if(!length(ID)==length(unique(ID))){
        fit <- tryCatch({
          fit1 <- lme(formula, data = dat_sub, random= ~ 1 | ID)
        }, error=function(err){
          fit1 <- try({lme(formula, data = dat_sub, random= ~ 1 | ID)}) 
          return(fit1)
        })
        
        ###################################
        # Summarize Coefficient Estimates #
        ###################################
        
        if (class(fit) != "try-error"){
          para<-as.data.frame(summary(fit)$tTable)[-1,-c(2:4)]
          colnames(para)<-c('coef', 'pval')
          para$metadata<-colnames(metadata)
          para$feature<-colnames(features)[x]
          rownames(para)<-NULL
        }
        else{
          print(paste("Fitting problem for feature", x, "returning NA"))
          para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
          colnames(para)<-c('coef', 'pval')
          para$metadata<-colnames(metadata)
          para$features<-colnames(features)[x]
          rownames(para)<-NULL
        } 
      }
    
    #######################
    # Fixed effects model #
    #######################
      
      else{
        fit <- tryCatch({
          fit1 <- glm(formula, data = dat_sub, family='gaussian')
        }, error=function(err){
          fit1 <- try({glm(formula, data = dat_sub, family='gaussian')}) 
          return(fit1)
        })
        
        ###################################
        # Summarize Coefficient Estimates #
        ###################################
        
        if (class(fit) != "try-error"){
          para<-as.data.frame(summary(fit)$coefficients)[-1,-c(2:3)]
          colnames(para)<-c('coef', 'pval')
          para$metadata<-colnames(metadata)
          para$feature<-colnames(features)[x]
          rownames(para)<-NULL
        } 
        else{
          print(paste("Fitting problem for feature", x, "returning NA"))
          para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
          colnames(para)<-c('coef', 'pval')
          para$metadata<-colnames(metadata)
          para$features<-colnames(features)[x]
          rownames(para)<-NULL
        }
      }
    
    return(para)
  })    
   
  ###################
  # Combine results #
  ###################
  
  paras<-do.call(rbind, paras)
  
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

################################
# Fit LM To A List of Datasets #
################################

list.LM<-function(physeq, transformation='NONE', multiple_qvalues = TRUE){
  foreach(physeq=physeq, .export=c("pvalueAdjustment_HM", "append_qvalues",
                                   "fit.LM", "AST", "LOG", "LOGIT"),
          .packages=c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                      "pbapply", "car", "nlme"), 
          .errorhandling = 'remove') %dopar% 
  {
    start.time <- Sys.time()
    features<-physeq$features; 
    metadata<-physeq$metadata;
    libSize <- physeq$libSize;
    ID<-physeq$ID;
    DD<-fit.LM(features, metadata, libSize, ID, transformation, multiple_qvalues)
    DD$pairwiseAssociation<-paste('pairwiseAssociation', 1:nrow(DD), sep='')
    wh.TP = intersect(grep("[[:print:]]+\\_TP$", DD$metadata), grep("[[:print:]]+\\_TP$", DD$feature))
    newname = paste0(DD$pairwiseAssociation[wh.TP], "_TP")
    DD$pairwiseAssociation[wh.TP] <- newname;
    DD<-dplyr::select(DD, c('pairwiseAssociation', 'feature', 'metadata'), everything())
    stop.time <- Sys.time()
    time<-as.numeric(round(difftime(stop.time, start.time, units="min"),3), units = "mins")
    DD$time<-time
    return(DD)
  }
}

