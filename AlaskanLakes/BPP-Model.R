# Compiled Files of Charcoal Process Point Model
# Aidan Draper

rm(list=ls())

# import libraries
library(spBayes)
library(ggplot2)
library(coda)
library(gridExtra)
library(mgcv)
library(zoo)
library(Matrix)
library(stringr)
library(rlecuyer)
library(cowplot)
# uncomment if you use Malcolm's snowfall approach
#library(snowfall)
#library(snow)


#STEP 1: Setup and Data Manipulation
setwd("~/Documents/Junior_Year/DISC_REU/DISC_bayesian_model/AlaskanLakes/")

char.dat = read.csv("chopperCHAR.csv",header=TRUE)
names(char.dat) = c("depth.a","depth.b","age.a","age.b","vol","count")

head(char.dat)

# create covariates to use in model                                 covariates:
char.dat$age = with(char.dat, (age.a+age.b)/2)                      # age of sediment core section
char.dat$sed.rate = with(char.dat, (depth.b-depth.a)/(age.b-age.a)) # sedimentation rate
char.dat$char = with(char.dat, (count/vol)*sed.rate)                # charcoal accumulation rate
char.dat$offset = with(char.dat,age.b-age.a)                        # time took for sediment core to form
char.dat$influx = with(char.dat,sed.rate/vol)                       # sediment influx rate?

# Drop last row of data with no age information
char.dat = char.dat[-nrow(char.dat),]

# plot 1: charcoal count over time
qplot(age,count,data=char.dat) + geom_smooth(method="loess", span=0.08) + theme_bw()


#STEP 2: Approximate background and foreground intensity (using cubic base spline) and set offset
char.dat$age.c = with(char.dat,age-min(age))
char.dat$age.s = with(char.dat,age.c/max(age.c))

n = nrow(char.dat)
n.knots = 101 # this variable will change for each model based on time interval

# note: smoothCon is part of the mgcv package for constructing the smooth terms in a GAM model
CRbasis = smoothCon(s(age.s,k=n.knots,bs="cr"),data=char.dat,knots=NULL,absorb.cons=TRUE,
                    scale.penalty=TRUE)
Sb = CRbasis[[1]]$S[[1]]
X = CRbasis[[1]]$X
knots = as.numeric(CRbasis[[1]]$xp)
S.scale = CRbasis[[1]]$S.scale

TT = char.dat$offset # sets offset (time interval between a and b calculated in STEP 1)


#STEP 3: Generate Starting Values (PRIOR?) and fit the first GAM models (m1 and m2) for Background and Foreground

# note: (rollapply is part of the zoo package) applies "quantile(x,0.9)" to get a rolling 90% threshold as more charcoal counts are included
# note: (na.approx is part of the zoo package) approximates the missing values

char.dat$t.hold = rollapply(char.dat$count,50,function(x)quantile(x,0.9),fill="extend")

# Separate charcoal counts
y.back = y.fore = char.dat$count
y.back[char.dat$count>char.dat$t.hold] = NA # removes count value for background if the value is distributed in the 90% quantile
y.fore[char.dat$count<char.dat$t.hold] = NA # removes count value for background if the value is NOT distributed in the 90% quantile

#approximate the missing values (background = na.approx, foreground = randomly pulls value from binomial distribution)
y.back = round(na.approx(y.back)) 
y.fore[is.na(y.fore)] = rbinom(sum(is.na(y.fore)),1,0.1) # ask about why we add a bunch of zeros instead of NAs

# Plot background and foreground counts
plot(char.dat$age,y.back)
plot(char.dat$age,y.fore)

# Add separate counts to data.frame
char.dat$count.b = y.back
char.dat$count.f = y.fore

# Fit GAM Model to background (gam is part of mgcv package)
m1 = gam(y.back~X,family=poisson,offset=log(TT),paraPen=list(X=list(Sb,sp=0.001*S.scale)))
# summary(m1)

# plots the points and our fitted GAM model for background
plot(char.dat$age,y.back,ylim=c(0,100))
lines(char.dat$age,fitted(m1),col="blue",lwd=2)

# Fit GAM Model to foreground
m2 = gam(y.fore~X,family=poisson,offset=log(TT),paraPen=list(X=list(Sb,sp=(1e-7)*S.scale)))
# summary(m2)

# plots the points and our fitted GAM model for foreground
plot(char.dat$age,y.fore)
lines(char.dat$age,fitted(m2),col="blue",lwd=2)


#STEP 4: Create your holdout set (25% of your data) - 
# uses gridded search to find both variances and decide holdout points

ho.idx = rep(NA,round(0.25*n)) # empty vector of 25 percent of data
ho.idx[1] = sample(1:n,1)
exclude = ho.idx[1] + c(-2,-1,0,1,2)

# decides what charcoal counts to exclude
for(i in 2:length(ho.idx)){
  check = TRUE
  while (check == TRUE){
    tmp = sample(1:n,1)
    check = tmp %in% exclude
  }
  ho.idx[i] = tmp
  exclude = c(exclude,tmp + c(-2,-1,0,1,2)) # will not exclude points +/- 2 of point 
}

# sorts the indexed hold out points
ho.idx = sort(ho.idx)

# plots charcoal count over time with red holdout points on top 
plot(char.dat$age,char.dat$count,pch=20,xlim=rev(range(char.dat$age)))
points(char.dat$age[ho.idx],char.dat$count[ho.idx],pch=20,col="red",cex=1.2)
legend(x=500,y=250,legend=c("Held Out Count"),pch=c(20),col=c("red"))

# distribution of the difference between two holdout points 
hist(diff(ho.idx))

# records the number of holdout points and new dataset length
n.hold = length(ho.idx)
n.obs = n - n.hold

# resets parameters
char.obs = char.dat[-ho.idx,] # creates new dataset for model
char.ho = char.dat[ho.idx,] # creates a dataset of our holdout points to check our model
y.obs = char.obs$count
y.oos = char.ho$count
X.obs = X[-ho.idx,]
X.ho = X[ho.idx,]
TT.obs = char.obs$offset
TT.ho = char.ho$offset


#STEP 5: Cross-validate your model against the holdout points to see the model's accuracy

# Define number of chains
q = 3

# Define background and foreground penalties (based on the smoothing scale of the splines)
n.pen = 5
sig2.b = seq(100,500,length.out=n.pen)/S.scale
sig2.f = seq(1e8,5e8,length.out=n.pen)/S.scale
pen.params = expand.grid(sig2.b=sig2.b,sig2.f=sig2.f)
G = nrow(pen.params)

# Define run index
run.idx = expand.grid(sig2.b=sig2.b,sig2.f=sig2.f,chain=1:q) # duplicate your (sigmoid function) penalties by the number of chains
run.idx = run.idx[,c("chain","sig2.b","sig2.f")]             # rearrange columns
run.idx$id = sapply(1:nrow(run.idx),function(i,x=run.idx)    # create your id variable
  paste(str_pad(which((1:q)%in%x$chain[i]),2,pad="0"),
        str_pad(which(sig2.b%in%x$sig2.b[i]),2,pad="0"),
        str_pad(which(sig2.f%in%x$sig2.f[i]),2,pad="0"),sep=""))



# FUNCTION (for later): Define log target density function
ld = function(theta,inputs,priors){
  
  for (i in 1:length(inputs)){
    tmp = names(inputs)
    assign(tmp[i],inputs[[i]])
  }
  
  for (i in 1:length(priors)){
    tmp = names(priors)
    assign(tmp[i],priors[[i]])
  }
  
  b0.b = theta[b0.b.idx]
  beta.b = theta[beta.b.idx]
  b0.f = theta[b0.f.idx]
  beta.f = theta[beta.f.idx]
  
  lam.Tb = exp(b0.b*one + X%*%beta.b)*TT
  lam.Tf = exp(b0.f*one + X%*%beta.f)*TT
  lam.T = lam.Tb + lam.Tf
  
  ld = sum(y*log(lam.T)) - sum(lam.T) -
    (1/(2*sig2.0))*crossprod(b0.b-mu.0) -
    (1/(2*sig2.b))*crossprod(beta.b-mu.b,S%*%(beta.b-mu.b)) -
    (1/(2*sig2.0))*crossprod(b0.f-mu.0) -
    (1/(2*sig2.f))*crossprod(beta.f-mu.f,S%*%(beta.f-mu.f))
  
  return(ld)
  
}



# Setup cross-validation inputs as a large list and save is to an .rda file

# create a list of the parameter names
param.names = c("b0.b",paste("beta.b[",seq(n.knots-1),"]",sep=""),"b0.f",
                paste("beta.f[",seq(n.knots-1),"]",sep=""))

# get index of column title
b0.b.idx = which(param.names=="b0.b")
b0.f.idx = which(param.names=="b0.f")

# get indecies where there is not a "column title"
beta.b.idx = which(param.names%in%paste("beta.b[",seq(n.knots-1),"]",sep=""))
beta.f.idx = which(param.names%in%paste("beta.f[",seq(n.knots-1),"]",sep=""))

mod.inputs = list() # creates new list

for (j in 1:q){ # for each chain
  for (g in 1:G){ # and for each penalty value 
    
    # idx = i + K*(g-1) + K*G*(j-1)
    
    idx = g + G*(j-1) # index is equal to the iterator (indexer) in the run.idx var - weird way to do this IMO
    
    # Define inputs - they all are the CRBasis (cubic spline results) from the oberservation model?
    y = y.obs
    TT = TT.obs
    one = rep(1,n.obs) # vector of 1s for n.obs (which is length(dataset) - length(holdout set))
    X = X.obs 
    inputs = list(y=y,S=Sb,one=one,X=X,TT=TT,b0.b.idx=b0.b.idx,beta.b.idx=beta.b.idx,
                  b0.f.idx=b0.f.idx,beta.f.idx=beta.f.idx)
    
    #### Define priors ####
    
    priors = list(mu.0=0,sig2.0=1e7,mu.b=rep(0,n.knots-1),mu.f=rep(0,n.knots-1),
                  sig2.b=pen.params[g,"sig2.b"],sig2.f=pen.params[g,"sig2.f"])
    
    # Define starting values
    b0.b = coef(m1)[1] + runif(1,-0.01,0.01)
    beta.b = coef(m1)[2:n.knots] + runif(n.knots-1,-0.01,0.01)
    b0.f = coef(m2)[1] + runif(1,-0.01,0.01)
    beta.f = coef(m2)[2:n.knots] + runif(n.knots-1,-0.01,0.01)
    
    
    starting = c(b0.b,beta.b,b0.f,beta.f)
    
    mod.inputs[[idx]] = list(ltd=ld,inputs=inputs,priors=priors,tune=0.1,
                             starting=starting)
    
  }
}

names(mod.inputs) = paste("run",run.idx$id,sep="")

#save(mod.inputs,file="InputsCVchopper.rda")
rm(mod.inputs) # to save space



# STEP 6: Actually run the cross-validation, works for multiple lakes

# Malcolm used Snowfall and ran the method within this file. However, we ran runCVbatch.R, runCVbatch.sh, and submitted a 
# job request to Notre Dame's HTCondor (cv_run.submit). If you do not have a grid system to use, you can try running the 
# code commented out below. Warning, it is computationally demanding and could run for ~24 or more depending on your
# computer, how many lakes and the size of the input.

# # FUNCTION: from the RunAdaptMCMC_CV.R file, runs the cross validation
# run.cv = function(x=mod.inputs){
#   
#   attach(x)
#   
#   out = try(adaptMetropGibbs(ltd=ltd,starting=starting,tuning=tune,batch=2400,report=400,
#                              inputs=inputs,priors=priors))
#   
#   n.sim = 2400*25
#   n.burn = 50000
#   n.step = 5
#   save.idx = seq(n.burn+1,n.sim,n.step)
#   
#   if (class(out)=="try-error"){
#     samps = array(NA,dim=c(length(save.idx),length(starting)))
#   } else {
#     samps = out$p.theta.samples[save.idx,] 
#   }
#   
#   return(samps)
#   
#   detach(x)
#   
# }

##### Original Code Tactic. However, we decided to use Notre Dames' HTCondor instead ####

# lakes = c("chopper")
# 
# for (i in 1:length(lakes)){
#   
#   lake.inputs = paste("InputsCV",lakes[i],".rda",sep="")
#   
#   
#   load(lake.inputs)
#   
#   #### THIS IS WHERE I LEFT OFF #### (need to figre out parallel computing to make it run faster)
#   out <- for(i in 1:idx) run.cv(mod.inputs[[i]]) # only way I could attach x properly so it could read variables
#   
#   #source("RunAdaptMCMC_CV.R")
# 
#   track.file = paste("CV",lakes[i],".txt",sep="")
# 
#   #sfInit(parallel=T,cpus=19,slaveOutfile=track.file)
#   #sfClusterSetupRNG()
# 
#   #out = sfClusterApplyLB(mod.inputs,run.cv)
# 
#   #sfStop()
# 
#   out.file = paste("OutCV",lakes[i],".rda",sep="")
#   
#   #save(out,file=out.file)
#   
# }

#### END OF ORIGINAL CROSS VALIDATION ####


# STEP 7: Process Cross Validation Output from HTCondor 
# (make sure to scp (remote copy) your output file to your personal computer first)


#### RUN THE BATCH FILE ON A CLUSTERED SYSTEM ####

# Then, Load CV results
load("cv_results_chopper.rda")
names(cv.results) = paste("run",run.idx$id,sep="")

chain.idx = t(sapply(1:G,function(i) i + G*(0:2)))

n.save = nrow(cv.results[[1]])

lam.T.sim = y.sim = list()

# Define MCMC objects for model parameters
for (g in 1:G){
  
  b0.b.sim = mcmc.list(lapply(1:q,function(j,x=cv.results,id=chain.idx[g,])
    mcmc(x[[id[j]]][,b0.b.idx])))
  
  beta.b.sim = mcmc.list(lapply(1:q,function(j,x=cv.results,id=chain.idx[g,])
    mcmc(x[[id[j]]][,beta.b.idx])))
  
  b0.f.sim = mcmc.list(lapply(1:q,function(j,x=cv.results,id=chain.idx[g,])
    mcmc(x[[id[j]]][,b0.f.idx])))
  
  beta.f.sim = mcmc.list(lapply(1:q,function(j,x=cv.results,id=chain.idx[g,])
    mcmc(x[[id[j]]][,beta.f.idx])))
  
  lam.Tb.sim = mcmc.list(lapply(1:q,function(j)
    mcmc(t(sapply(1:n.save,function(i,b0=b0.b.sim[[j]],b=beta.b.sim[[j]],x=X.ho,off=TT.ho)
      exp(b0[i] + x%*%b[i,])*off)))))
  
  lam.Tf.sim = mcmc.list(lapply(1:q,function(j)
    mcmc(t(sapply(1:n.save,function(i,b0=b0.f.sim[[j]],b=beta.f.sim[[j]],x=X.ho,off=TT.ho)
      exp(b0[i] + x%*%b[i,])*off)))))
  
  lam.T.sim[[g]] = mcmc.list(lapply(1:q,function(j)
    mcmc(t(sapply(1:n.save,function(i,lb=lam.Tb.sim[[j]],lf=lam.Tf.sim[[j]])
      lb[i,]+lf[i,])))))
  
  y.sim[[g]] = mcmc.list(lapply(1:q,function(j)
    mcmc(t(sapply(1:n.save,function(i,lam=lam.T.sim[[g]][[j]],n=n.hold)
      rpois(n,lam[i,]))))))
  
}

# Derive posterior means and variances for each held out charcoal count
post.mean.y = t(sapply(y.sim,function(x)apply(do.call(rbind,x),2,mean)))
post.var.y = t(sapply(y.sim,function(x)apply(do.call(rbind,x),2,var)))

# Plot posterior mean counts for holdout set (commented out here for simplicity)

#x11()
#par(ask=T)
#for (i in 1:G){
#plot(char.dat$age,char.dat$count,"l")
#points(char.dat$age[ho.idx],post.mean.y[i,],pch=20,col="red",cex=1.2)
#}

post.summ = run.idx[1:G,c("id","sig2.b","sig2.f")]

# Estimate posterior predictive loss ("ppl") for each set of penalty values
post.summ$ppl = sapply(1:G,function(i,yo=y.oos,ymean=post.mean.y,yvar=post.var.y)
  crossprod(yo-ymean[i,]) + sum(yvar[i,]))

# Order the penalty sets in ascending order of ppl (lower values are preferred, 
# so the first row entry of the post.summ array indicates the optimal penalty values)
post.summ = post.summ[order(post.summ$ppl),]
post.summ


# STEP 8: Setup inputs for final run of model using optimal penalties

# Model inputs
X = CRbasis[[1]]$X # this may be incorrect
y = char.dat$count
TT = char.dat$offset
age = char.dat$age
one = rep(1,n)

inputs = list(y=y,S=Sb,one=one,X=X,TT=TT,b0.b.idx=b0.b.idx,beta.b.idx=beta.b.idx,
              b0.f.idx=b0.f.idx,beta.f.idx=beta.f.idx)

sig2.b = post.summ[1,"sig2.b"]
sig2.f = post.summ[1,"sig2.f"]

priors = list(mu.0=0,sig2.0=1e7,mu.b=rep(0,n.knots-1),mu.f=rep(0,n.knots-1),
              sig2.b=sig2.b,sig2.f=sig2.f)

q = 3

mod.inputs = list()

#set.seed(24678)
for (i in 1:q){
  
  # Define starting values
  b0.b = coef(m1)[1] + runif(1,-0.01,0.01)
  beta.b = coef(m1)[2:n.knots] + runif(n.knots-1,-0.01,0.01)
  b0.f = coef(m2)[1] + runif(1,-0.01,0.01)
  beta.f = coef(m2)[2:n.knots] + runif(n.knots-1,-0.01,0.01)
  
  starting = c(b0.b,beta.b,b0.f,beta.f)
  
  mod.inputs[[i]] = list(ltd=ld,inputs=inputs,priors=priors,tune=0.1,
                         starting=starting)
  
}

#save(mod.inputs,file="InputsMCRchopper.rda")




###### Step 9: PROCESS FINAL MODEL ######
load("~/Documents/Junior_Year/DISC_REU/DISC_bayesian_model/AlaskanLakes/mcmc_results_chopper.rda")

q = length(mcmc.results)

n.samp = nrow(mcmc.results[[1]]$p.theta.samples)
n.burn = n.samp * .25
n.step = 10
save.idx = seq((n.burn+1): n.samp) + n.burn # weird glitch in seq() function...
save.idx <- save.idx[which(save.idx %% 10 == 1)]
n.save = length(save.idx)

b0.b.sim = mcmc.list(lapply(mcmc.results,function(x)mcmc(x$p.theta.samples[save.idx,b0.b.idx])))
beta.b.sim = mcmc.list(lapply(mcmc.results,function(x)mcmc(x$p.theta.samples[save.idx,beta.b.idx])))
b0.f.sim = mcmc.list(lapply(mcmc.results,function(x)mcmc(x$p.theta.samples[save.idx,b0.f.idx])))
beta.f.sim = mcmc.list(lapply(mcmc.results,function(x)mcmc(x$p.theta.samples[save.idx,beta.f.idx])))

X = CRbasis[[1]]$X

lam.Tb.sim = mcmc.list(lapply(1:q,function(j)
  mcmc(t(sapply(1:n.save,function(i,b0=b0.b.sim[[j]],b=beta.b.sim[[j]],x=X,off=TT)
    exp(b0[i] + x%*%b[i,])*off)))))

lam.Tf.sim = mcmc.list(lapply(1:q,function(j)
  mcmc(t(sapply(1:n.save,function(i,b0=b0.f.sim[[j]],b=beta.f.sim[[j]],x=X,off=TT)
    exp(b0[i] + x%*%b[i,])*off)))))

lam.T.sim = mcmc.list(lapply(1:q,function(j)
  mcmc(t(sapply(1:n.save,function(i,lb=lam.Tb.sim[[j]],lf=lam.Tf.sim[[j]])
    lb[i,]+lf[i,])))))

pT.sim = mcmc.list(lapply(1:q,function(j)
  mcmc(t(sapply(1:n.save,function(i,lb=lam.Tb.sim[[j]],lf=lam.Tf.sim[[j]])
    lf[i,]/(lb[i,]+lf[i,]))))))

# Check convergence of probability of fire parameters
# (commented mcmc.results for simplicity)

#x11(width=5,height=8)
#plot(pT.sim,ask=T)
# 
#gelman.diag(pT.sim)

# Fitted charcoal counts subplot
summ.lam.T = summary(lam.T.sim)
char.dat$mu.post.mean = summ.lam.T$statistics[,"Mean"]
p1 = ggplot(aes(x=age,y=count), data=char.dat) +
  geom_point(size=1.5,shape=1) +
  geom_line(aes(x=age,y=mu.post.mean,color="Fitted Charcoal Count"),size=1) +
  scale_color_manual(values=c("blue")) +
  ylab("Count") + scale_x_reverse() + xlab("") +
  theme_bw() + theme(legend.title=element_blank(),
                     legend.position=c(0.15,0.8),
                     panel.grid=element_blank())

# Background/foreground intensity subplot
lam.b.sim = lam.f.sim = array(NA,dim=c(q*n.save,n))

for (j in 1:q){
  
  lam.b.sim[(1:n.save)+n.save*(j-1),] = 
    t(sapply(1:n.save,function(i,b0=b0.b.sim[[j]],x=X,B=beta.b.sim[[j]])
      exp(b0[i] + x%*%B[i,])))
  
  lam.f.sim[(1:n.save)+n.save*(j-1),] = 
    t(sapply(1:n.save,function(i,b0=b0.f.sim[[j]],x=X,B=beta.f.sim[[j]])
      exp(b0[i] + x%*%B[i,])))
  
}

char.dat$lam.b.post.mean = apply(lam.b.sim,2,mean)
char.dat$lam.f.post.mean = apply(lam.f.sim,2,mean)
p2 = ggplot(aes(x=age,y=lam.b.post.mean), data=char.dat) +
  geom_line(aes(color="background"),size=1,alpha=I(0.75)) +
  geom_line(aes(x=age,y=lam.f.post.mean,color="fire"),size=1,alpha=I(0.75)) +
  scale_color_manual(values=c("green","red")) +
  ylab("Intensity") + scale_x_reverse() + xlab("") +
  theme_bw() + theme(legend.title=element_blank(),
                     legend.position=c(0.1,0.7),
                     panel.grid=element_blank())

# Probability of fire subplot
summ.pT = summary(pT.sim)
char.dat$pf = summ.pT$statistics[,"Mean"]
char.dat$pf.lower = summ.pT$quantiles[,"2.5%"]
char.dat$pf.upper = summ.pT$quantiles[,"97.5%"]
p4 = ggplot(aes(x=age,y=pf),data=char.dat) +
  geom_ribbon(aes(x=age,ymin=pf.lower,ymax=pf.upper),fill="grey") +
  geom_line(size=0.8) +
  # geom_line(aes(y=pT),color="purple") +
  #   geom_point() +
  geom_hline(aes(yintercept=0.9),color="red") +
  ylab("Prob(Fire)") + scale_x_reverse() + xlab("Time (YBP)") +
  theme_bw() +
  theme(panel.grid=element_blank())


print(plot_grid(p1,p2,p4,nrow=3,align="v"))

