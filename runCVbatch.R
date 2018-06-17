# Aidan Draper
# June 12, 2018
# Batch system file for cross-validation of charcoal Bayesian model

install.packages('spBayes')
library(spBayes)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("Three arguments must be supplied (lake name, input process). Run again.", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  lake.name = args[1]
  process = int(args[2]) - 1 
  if (process==0){
    start = process + 1
    stop = process + 3
  } else {
    start = process * 3
    stop = (process * 3) + 2
  }
} else {
  stop("Too many arguments supplied. Run again.", call.=FALSE)
}

lake.inputs = paste("InputsCV", lake.name, ".rda",sep="")

load(lake.inputs)


###############################################################################
## Function created by Malcolm Itter
## May 25, 2018
###############################################################################
## RunAdaptMCMC_CV - Function to run MCMC simulation for univariate model
## using cross validation inputs
###############################################################################

run.cv = function(x){
  
  attach(x)
  
  out = try(adaptMetropGibbs(ltd=ltd,starting=starting,tuning=tune,batch=2400,report=400,
                             inputs=inputs,priors=priors))
  
  n.sim = 2400*25
  n.burn = 50000
  n.step = 5
  save.idx = seq(n.burn+1,n.sim,n.step)
  
  if (class(out)=="try-error"){
    samps = array(NA,dim=c(length(save.idx),length(starting)))
  } else {
    samps = out$p.theta.samples[save.idx,] 
  }
  
  return(samps)
  
  detach(x)
  
}

####### END OF Malcolm's CODE ######


results = vector("list", ((stop - start) + 1))

for (i in start+1 : stop){
  results[[i]] = run.cv(mod.inputs[i])
}

return(results)

