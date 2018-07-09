# Aidan Draper
# June 29,2018

# this file takes the output .rda files from the HTCondor system process and produces one final .rda file called 
# cv_results_{lake name}.rda

#change the directory to your locally copied output file from condor
setwd('~/Documents/Junior_Year/DISC_REU/DISC_bayesian_model/BigWoodsLakes/outputBW')

# make sure to adjust your parameters to match your HTCondor job output
proc <- 25 # number of processes
inp <- 3 # number of CVs run per process
lake.name <- "crystal"

cv.results = vector("list", ((proc) * inp))

count <- 1
for(x in 0 : (proc - 1)){
  load(paste('OutCV',lake.name,x,'.rda',sep=''))
  #paste(results, '\n')
  posit <- count
  for (y in 1 : length(results)){
    if (length(results[[y]]) != 0){
      cv.results[[count]] <- results[[y]]
      count <- count + 1
    }
  }
  rm(results) # to make sure we are grabbing the new results input
}
setwd('~/Documents/Junior_Year/DISC_REU/DISC_bayesian_model/BigWoodsLakes/')
save(cv.results,file=paste('cv_results_',lake.name,'.rda', sep=""))
