#!/usr/bin/env Rscript

# Install and load Packages
if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org') 
suppressPackageStartupMessages(library("pacman"))
if(! require("devtools")) install.packages("devtools", repos='http://cran.us.r-project.org') 
suppressPackageStartupMessages(library("devtools"))
if(! require("sparseDOSSA")) devtools::install_github('biobakery/sparseDOSSA@varyLibSize')
suppressPackageStartupMessages(library("sparseDOSSA"))
pacman::p_load('pkgmaker', 'optparse', 'parallel', 'stringi', 'MASS', 'doParallel')
suppressPackageStartupMessages(library("pkgmaker"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("stringi"))
suppressPackageStartupMessages(library("doParallel"))
 
# Command Line Usage
option_list = list(
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
  c("-q", "--readDepth"), default=50000, # q stands for quantity and quality of sequencing depth or library size
  type = "integer"),
make_option(
  c("-g", "--noParallel"), default=FALSE, # g stands for grid
  action = "store_true"),
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
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)
print(opt)

# Set Parameters for Testing Purposes
# noZeroInflate<- FALSE # High-level parameter
# RandomEffect<-FALSE # High-level parameter
# metadataType <- 'MVA' # High-level parameter
# nSubjects <- 50 # Low-level parameter
# nPerSubject<-1 # Low-level parameter
# nMicrobes <- 100 # Low-level parameter
# spikeMicrobes <- 0.1 # Low-level parameter
# nMetadata<- 5 # Low-level parameter
# spikeMetadata<- 0.2 # Low-level parameter
# effectSize<- 5 # Low-level parameter
# readDepth<- 50000 # Default parameter
# noParallel<- FALSE # Default parameter
# nIterations<- 100 # Default parameter
# rSeed<- 1234 # Default parameter
# nCores<- 4 # Default parameter
# workingDirectory <- '/Users/hmallick/Dropbox (Huttenhower Lab)/Maaslin2' # Default parameter

# Extract Parameters
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
noParallel<-opt$options$noParallel # Default parameter
nIterations<- opt$options$nIterations # Default parameter
rSeed<- opt$options$rSeed # Default parameter
nCores<- opt$options$nCores # Default parameter
workingDirectory <- opt$options$workingDirectory # Default parameter

# Load Utility Functions
pkgmaker::source_files(paste(workingDirectory, '/Library/sparseDOSSAhelperFunctions.R', sep=''))


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


# If input does not exist, then only generate data. Otherwise, skip.

if (!file.exists(inputStringR)){
  # Do some operation based on user input
  simlist<-trigger_sparseDOSSA_Simulator(noZeroInflate=noZeroInflate,
                                         RandomEffect=RandomEffect,
                                         metadataType=metadataType,
                                         nSubjects=nSubjects,
                                         nPerSubject=nPerSubject,
                                         nMicrobes=nMicrobes,
                                         spikeMicrobes = spikeMicrobes, 
                                         nMetadata = nMetadata, 
                                         spikeMetadata=spikeMetadata,
                                         effectSize = effectSize,
                                         readDepth = readDepth,
                                         noParallel = noParallel,
                                         nIterations=nIterations,
                                         rSeed=rSeed,
                                         nCores=nCores)
  
  
  
  
  save(simlist, file=inputStringR)
} else{
  print("Input file already exists. No new data generated!")
  load(inputStringR)
}     


#######################################
# Delete Temporary sparseDOSSA Files  #
#######################################

if (file.exists("SyntheticMicrobiome-Counts.pcl")) file.remove("SyntheticMicrobiome-Counts.pcl")
if (file.exists("SyntheticMicrobiome.pcl")) file.remove("SyntheticMicrobiome.pcl")
if (file.exists("SyntheticMicrobiomeParameterFile.txt")) file.remove("SyntheticMicrobiomeParameterFile.txt")
