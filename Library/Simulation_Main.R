#!/usr/bin/env Rscript

# Install and Load Packages
if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org') 
suppressPackageStartupMessages(library("pacman"))
if(! require("devtools")) install.packages("devtools", repos='http://cran.us.r-project.org') 
suppressPackageStartupMessages(library("devtools"))
pacman::p_load('pkgmaker', 'optparse', 'parallel', 'stringi', 'vegan', 'doParallel', 'ROCR', 'plyr', 'chemometrics')
if(!require('phyloseq')){
  source('http://bioconductor.org/biocLite.R')
  biocLite('phyloseq')
}
suppressPackageStartupMessages(library("pkgmaker"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("stringi"))
suppressPackageStartupMessages(library("vegan"))
suppressPackageStartupMessages(library("doParallel"))
suppressPackageStartupMessages(library("ROCR"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("edgeR")) # TMM
suppressPackageStartupMessages(library("metagenomeSeq")) # CSS
suppressPackageStartupMessages(library("DESeq2")) # RLE
suppressPackageStartupMessages(library("phyloseq")) # RFA


# Command Line Usage
option_list = list(
  make_option(
    c("-a", "--methodName"), # a stands for association method
    type = "character"),
  make_option(
    c("-z", "--noZeroInflate"), default=FALSE, # z stands for zero-inflation
    action = "store_true"),
  make_option(
    c("-l", "--RandomEffect"), default=FALSE, # l stands for longitudinal design
    action = "store_true"),
  make_option(
    c("-d", "--metadataType"), # d stands for design matrix
    type = "character"),
  make_option(
    c("-n", "--nSubjects"), # n stands for number of subjects
    type = "integer"),
  make_option(
    c("-b", "--nPerSubject"), default=1), # b stands for block size
  make_option(
    c("-f", "--nMicrobes"), # f stands for number of features
    type = "integer"),
  make_option(
    c("-s", "--spikeMicrobes"), # s stands for how sparsely the features are associated
    type = "numeric"),
  make_option(
    c("-m", "--nMetadata"), # m stands for number of metadata
    type = "integer"),
  make_option(
    c("-p", "--spikeMetadata"), # p stands for what percentage of weights are non-zero
    type = "numeric"),
  make_option(
    c("-e", "--effectSize"), # e stands for effect size
    type = "integer"),
  make_option(
    c("-q", "--readDepth"), default=50000, # q stands for quantity of sequencing depth or library size
    type = "integer"),
  make_option(
    c("-t", "--nIterations"), default=100, # t stands for how many times an experiment is repeated
    type = "integer"),
  make_option(
    c("-r", "--rSeed"), default=1234, # r stands for reproducibility index (random seed)
    type = "integer"),
  make_option(
    c("-c", "--nCores"), default=4, # c stands for how many cores to be used in the analysis
    type = "integer"),
  make_option(
    c("-w", "--workingDirectory"), # w stands for working directory
    type = "character")) 

# Print Progress Message
cat("Running the following combination of method and parameters:", "\n");
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)
print(opt)

# Extract Parameters
methodName <- opt$options$methodName # High-level parameter
noZeroInflate<- opt$options$noZeroInflate # High-level parameter
RandomEffect<-opt$options$RandomEffect # High-level parameter
metadataType <- opt$options$metadataType # High-level parameter
nSubjects <- opt$options$nSubjects # Low-level parameter
nPerSubject<-opt$options$nPerSubject # Low-level parameter
nMicrobes <- opt$options$nMicrobes # Low-level parameter
spikeMicrobes <- opt$options$spikeMicrobes # Low-level parameter
nMetadata<- opt$options$nMetadata # Low-level parameter
spikeMetadata<- opt$options$spikeMetadata # Low-level parameter
effectSize<- opt$options$effectSize # Low-level parameter
readDepth<- opt$options$readDepth # Default parameter
nIterations<- opt$options$nIterations # Default parameter
rSeed<- opt$options$rSeed # Default parameter
nCores<- opt$options$nCores # Default parameter
workingDirectory <- opt$options$workingDirectory # Default parameter

# Set Parameters for Testing Purposes
# methodName = 'LM.CLR' # High-level parameter
# noZeroInflate<- FALSE # High-level parameter
# RandomEffect<-FALSE # High-level parameter
# metadataType <- 'UVB' # High-level parameter
# nSubjects <- 200 # Low-level parameter
# nPerSubject<-1 # Low-level parameter
# nMicrobes <- 500 # Low-level parameter
# spikeMicrobes <- 0.1 # Low-level parameter
# nMetadata<- 1 # Low-level parameter
# spikeMetadata<- 1 # Low-level parameter
# effectSize<- 10 # Low-level parameter
# readDepth<- 50000 # Default parameter
# nIterations<- 10 # Default parameter
# rSeed<- 1234 # Default parameter
# nCores<- 4 # Default parameter
# workingDirectory <- '/Users/hmallick/Dropbox (Huttenhower Lab)/Maaslin2' # Default parameter

# Load Utility Functions
pkgmaker::source_files(paste(workingDirectory, '/Library', sep=''), '*.R')

###########################################################
# ALL AVAILABLE MODELS + NORMALIZATIONS + TRANSFORMATIONS #
###########################################################

modelParams<-c('ANCOM', 'CPLM', 'DESeq2', 'edgeR', 'limma', 'limma2', 'limmaVOOM',
               'LM', 'LM2', 'metagenomeSeq','metagenomeSeq2', 'negbin', 'Spearman', 
               'Wilcoxon', 'ZIB','ZICP', 'ZINB')
normParams<-c('CSS', 'CLR', 'RLE', 'TMM', 'TSS')
transfParams<-c('AST', 'LOG', 'LOGIT')

###################################################
# Extract Method + Normalization + Transformation #
###################################################

inputMethodParams<-unlist(strsplit(methodName, "[.]"))
modelName<-inputMethodParams[inputMethodParams %in% modelParams]
normMethod<- inputMethodParams[inputMethodParams %in% normParams]
transfMethod<- inputMethodParams[inputMethodParams %in% transfParams]

if(is.character(modelName) & length(modelName) == 0) stop('The input model is invalid!')
if(is.character(transfMethod) & length(transfMethod) == 0) transfMethod='NONE'
if(is.character(normMethod) & length(normMethod) == 0) normMethod='NONE'


# Track Start Time
cat(c("Job started at:",date()), "\n")
start.time <- Sys.time()

# Extract Pre-generated Dataset 
if (noZeroInflate==TRUE && RandomEffect==TRUE){
  inputSubString<-'/noZeroInflate_RandomEffect'
}

if (noZeroInflate==TRUE && RandomEffect==FALSE){
  inputSubString<-'/noZeroInflate_noRandomEffect'
}

if (noZeroInflate==FALSE && RandomEffect==TRUE){
  inputSubString<-'/ZeroInflate_RandomEffect'
}

if (noZeroInflate==FALSE && RandomEffect==FALSE){
  inputSubString<-'/ZeroInflate_noRandomEffect'
}

inputString<-paste(inputSubString, metadataType, nSubjects, nPerSubject, nMicrobes, spikeMicrobes, nMetadata, spikeMetadata, effectSize, readDepth, sep='_')

# Create Input Directory
inputDirectory <- file.path(workingDirectory, 'Input')

if (!dir.exists(inputDirectory)){
  print("Input directory created!")
  dir.create(inputDirectory)
} else {
  print("Input directory already exists!")
}

inputStringR<-paste(inputDirectory, paste(inputString, '.RData', sep=''), sep='')

# Load Dataset
if (!file.exists(inputStringR)) stop('The input file does not exist. Generate the dataset first.') 
load(inputStringR)

# Set Up Clustering Environment
no_cores <- nCores  # detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

####################
# Simple Filtering #
####################

Threshold_Abundance = 0
Threshold_Prevalence = 0.1
simlist.filtered <- list.SimpleFilter(simlist, Threshold_Abundance, Threshold_Prevalence)

# Choose Appropriate Normalized Dataset 
if (normMethod=='NONE') {
  fitlist<-simlist.filtered
}

#####################
# RFA Normalization #
#####################

if (normMethod=='RFA') {
  RFAlist<-list.RFA(simlist.filtered); names(RFAlist) <- names(simlist.filtered)
  fitlist<-RFAlist
}

#####################
# TSS Normalization #
#####################

if (normMethod=='TSS') {
  TSSlist<-list.TSS(simlist.filtered); names(TSSlist) <- names(simlist.filtered)
  fitlist<-TSSlist
}

#####################
# CLR Normalization #
#####################
if (normMethod=='CLR') {
  CLRlist<-list.CLR(simlist.filtered); names(CLRlist) <- names(simlist.filtered)
  fitlist<-CLRlist
}

#####################
# CSS Normalization #
#####################

if (normMethod=='CSS') {
  CSSlist<-list.CSS(simlist.filtered); names(CSSlist) <-  names(simlist.filtered)
  fitlist<-CSSlist
}

#####################
# TMM Normalization #
#####################

if (normMethod=='TMM') {
  TMMlist<-list.TMM(simlist.filtered); names(TMMlist) <-  names(simlist.filtered)
  fitlist<-TMMlist
}

#####################
# RLE Normalization #
#####################

if (normMethod=='RLE') {
  RLElist<-list.RLE(simlist.filtered); names(simlist.filtered) <-  names(simlist.filtered)
  fitlist<-RLElist
}


# How Many Datasets to Run In A List 
reps<-1:nIterations
fitlist<-fitlist[reps]

# Assign Consistent Names to Results
if (noZeroInflate==TRUE && RandomEffect==TRUE){
  outputSubString<-'noZeroInflate_RandomEffect'
}

if (noZeroInflate==TRUE && RandomEffect==FALSE){
  outputSubString<-'noZeroInflate_noRandomEffect'
}

if (noZeroInflate==FALSE && RandomEffect==TRUE){
  outputSubString<-'ZeroInflate_RandomEffect'
}

if (noZeroInflate==FALSE && RandomEffect==FALSE){
  outputSubString<-'ZeroInflate_noRandomEffect'
}

outputString<-paste(methodName, outputSubString, metadataType, nSubjects, nPerSubject, nMicrobes, spikeMicrobes, nMetadata, spikeMetadata, effectSize, readDepth, sep='_')

# Create Output Directory
outputDirectory <- file.path(workingDirectory, 'Output')

if (!dir.exists(outputDirectory)){
  print("Output directory created!")
  dir.create(outputDirectory)
} else {
  print("Output directory already exists!")
}

# Save the Results in A .RData File
outputStringR<-paste(outputDirectory, paste(outputString, '.RData', sep=''), sep='/')

# Run the Model Only if the Output Does Not Exist
if (!file.exists(outputStringR)){
  
  # Model
  list.function<-eval(parse(text=paste('list', modelName, sep ='.')))
  output<-eval(parse(text="list.function(fitlist, transformation=transfMethod)")) 
 
  # Discard Failed Methods 
  if (length(output)<1) stop('Consistent error in the model fitting. No output returned.')
  names(output)<-paste(outputString, 1:length(output), sep='_')
  

  # Organize into A Coherent Data Frame
  dflist <- lapply(output, eval_res_list)
  names(dflist)<- names(output)
  simparamslabels = c("methodName", "ZeroInflate", "RandomEffect","metadataType", 
                      "nSubjects", "nPerSubject", "nMicrobes", "spikeMicrobes", "nMetadata", 
                      "spikeMetadata", "effectSize", "readDepth", "rep")
  
  df = make_power_df(dflist, simparamslabels)
  save(df, file=outputStringR)
} else{
  print("Output file already exists. No new results generated!")
  load(outputStringR)
}

# Stop the Cluster 
stopCluster(cl)


# Track End Time
stop.time <- Sys.time()
time<-round(difftime(stop.time, start.time, units="min"),3)
cat(c("Job finished at:",date()), "\n");
cat("Computational time:",time,"minutes \n")
cat("The output is in:", outputDirectory, fill=TRUE)
