#!/usr/bin/env Rscript
# Aidan Draper
# June 12, 2018
# Batch system file for MCMC simulation of charcoal Bayesian model

library(spBayes, lib = "/afs/crc.nd.edu/user/a/adraper2/condor/Rlibs")
library(mgcv, lib = "/afs/crc.nd.edu/user/a/adraper2/condor/Rlibs")

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Two arguments must be supplied (lake name, chains). Run again.", call.=FALSE)
} 

lake.name = args[1]
index <- as.integer(args[2]) + 1

lake.inputs = paste("InputsMCR", lake.name, ".rda",sep="")

load(lake.inputs)


###############################################################################
## Function created by Malcolm Itter
## May 25, 2018
###############################################################################
## RunAdaptMCMC_CV - Function to run MCMC simulation for univariate model
## using cross validation inputs
###############################################################################

run.adapt.mcmc = function(x){
  
  
  attach(x)
  
  out = adaptMetropGibbs(ltd=ltd,starting=starting,tuning=tune,batch=400,report=100,
                         inputs=inputs,priors=priors)
  
  detach(x)
  
  return(out)
  
}

####### END OF Malcolm's CODE ######

results = vector("list", 1)

results[[i]] = run.adapt.mcmc(mod.inputs[[index]])


out.file = paste("OutMCR",lake.name,process,".rda",sep="")

save(results,file=out.file)

