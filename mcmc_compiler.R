# Aidan Draper
# June 29,2018

# this file takes the output .rda files from the HTCondor system process and produces one final .rda file called 
# cv_results_{lake name}.rda

#change the directory to your locally copied output file from condor
setwd('~/Documents/Junior_Year/DISC_REU/DISC_bayesian_model/BigWoodsLakes/MCRrda')

# make sure to adjust your parameters to match your HTCondor job output
proc <- 3 # number of processes
lake.name <- "crystal"

mcmc.results = vector("list", proc)

for(x in 1 : proc){
  load(paste('OutMCR', lake.name, x,'.rda',sep=''))
  #paste(results, '\n')
  mcmc.results[[x]] <- results[[1]]
  rm(results) # to make sure we are grabbing the new results input
}

setwd('~/Documents/Junior_Year/DISC_REU/DISC_bayesian_model/BigWoodsLakes/')
save(mcmc.results,file=paste('mcmc_results_',lake.name,'.rda', sep=""))
