############################
# Zero-inflated Beta (ZIB) #
############################

###########################
# Load Essential Packages #
###########################

# pacman, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
load_essential_packages()

###########################################
# Load Dedicated Method-specific Packages #
###########################################

pacman::p_load('gamlss', 'pbapply', 'gamlss.dist', 'vegan')
if(! require("ZIBR")) {
  library(devtools)
  devtools::install_github("chvlyl/ZIBR")
}
library(ZIBR)

########################
# Fit ZIB To A Dataset #
########################

fit.ZIB <- function(features, 
                    metadata, 
                    libSize, 
                    ID, 
                    transformation,
                    multiple_qvalues){
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default ZIB model. Use NONE.')
  
  # TSS Normalizing the Data (ZIB assumes TSS-normalized or proportional data)
  features <- decostand(features, method="total", MARGIN=1)
  
  paras <- pbapply::pbsapply(1:ncol(features), simplify=FALSE, function(x){
  
    
    featuresVector <- features[, x]
    
    # Scrap All-Zero Features
    if(sum(featuresVector!=0)<1){
      print(paste("Cannot fit model to all zeroes for feature", x, "returning NA"))
      para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
      colnames(para)<-c('coef', 'pval')
      para$metadata<-colnames(metadata)
      para$features<-colnames(features)[x]
      rownames(para)<-NULL
    }  
    
    else{
    
      # Random effect adjustment
      if(!length(ID)==length(unique(ID))){
        nPerSubject<-round(length(ID)/length(unique(ID)))
        time.ind<-rep(1:nPerSubject, length(unique(ID)))

        # Fit Model
        fit <- tryCatch({
          fit1 <- ZIBR::zibr(logistic.cov=metadata, beta.cov=metadata, 
                             Y = featuresVector, subject.ind=ID, time.ind=time.ind)
        }, error=function(err){
          fit1 <- try({ZIBR::zibr(logistic.cov=metadata, beta.cov=metadata, 
                                  Y = featuresVector, subject.ind=ID, time.ind=time.ind)}) 
          return(fit1)
        })
        
        if (class(fit) != "try-error"){
          para<-as.data.frame(fit$beta.est.table)[-1,]
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
          para$feature<-colnames(features)[x]
          rownames(para)<-NULL
        }
      } else{
        
        # Fit Model
        dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata, libSize)
        formula<-as.formula(paste("expr ~ ", paste(colnames(metadata), collapse= "+")))
        fit <- tryCatch({
          fit1 <- gamlss(formula, data=dat_sub, family=BEZI(sigma.link="identity"), sigma.fo=~(libSize-1),trace=FALSE,control = gamlss.control(n.cyc = 100))
        }, error=function(err){
          fit1 <- try({gamlss(formula, data=dat_sub, family=BEZI(sigma.link="identity"), sigma.fo=~(libSize-1),trace=FALSE,control = gamlss.control(n.cyc = 100))}) 
          return(fit1)
        })
        
        if (class(fit) != "try-error"){
          para<-as.data.frame(summary(fit))[-c(1, ncol(metadata) + (2:3)),-c(2:3)]
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
          para$feature<-colnames(features)[x]
          rownames(para)<-NULL
        }
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

#################################
# Fit ZIB To A List of Datasets #
#################################

list.ZIB<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.ZIB"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        'gamlss', 'pbapply', 'gamlss.dist', 'vegan', 'ZIBR'),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.ZIB(features, metadata, libSize, ID, transformation, multiple_qvalues)
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

