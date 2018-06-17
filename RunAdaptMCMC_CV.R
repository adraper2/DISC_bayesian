###############################################################################
## Malcolm Itter
## May 25, 2018
###############################################################################
## RunAdaptMCMC_CV - Function to run MCMC simulation for univariate model
## using cross validation inputs
###############################################################################

run.cv = function(x){
  
  library(spBayes)
  
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
