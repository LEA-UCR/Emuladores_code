# if you are working local, setwd in simulation2 first!
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  i<-1
  type<-"Matern"
  model<-"SVI"
  analysis<-"M3"
  datasetfile=paste0("sim_data/dataset",model,type,i,".Rdata")
} else {
  i<-args[1]
  type<-args[2]
  model<-args[3]
  analysis<-args[4]
  datasetfile=paste0("sim_data/dataset",model,type,i,".Rdata")
}

# i<-1:100
# type<-'Exponential', "Matern"
# model<-'SVC', "SVI"
# analysis<-"M1: likelihood", "M2: Banerjee", "M3: FSA"

source("1.MRA_resolution_general.R")
source('covariances.R')
source('likelihoodK_general.R')
library(MCMCpack)
library(truncdist)
library(invgamma)

aa<-gen_resolution(datasetfile)
bordes<-aa[[1]];indicesW<-aa[[2]];knotsMRA<-aa[[3]]
nn<-aa[[4]];hh<-aa[[5]]

#### Metropolis Hastings

# initial values

phi <- 0.9
beta0 <- 0
beta1 <- 2
nu <- 1

startvalue <- c(phi,beta0,beta1,nu)
N <- dim(hh)[1]
npar <- length(startvalue)

# fixed
taub <- 1
taue <- 5

##################
# L functions    #
##################

f <- function(param) {
  phi   <- param[1]
  beta0 <- param[2]
  beta1 <- param[3]
  nu    <- param[4]
  sigma2 <- 1/taub
  if (analysis=="M1"){
    loglike <- likelihoodGaussian(nu,phi,beta0,
                  beta1,sigma2,taue,model,type)
    }else{
      if (analysis=="M2"){
    loglike <- likelihoodBanerjee(nu,phi,beta0,
                  beta1,sigma2,taue,model,type)
      }else {
    loglike <- likelihoodFSA_Block(nu,phi,beta0,
                  beta1,sigma2,taue,model,type)}}
  #logpriortaue <- (dgamma(taub,shape=0.5, scale=2, log=T))
  #logpriortaub <- dgamma(taub,shape=5, scale=2, log=T)
  logpriorphi   <- dunif(phi,0.80,1.00,log=TRUE) 
  logpriorbeta0 <- dnorm(0,1,log=TRUE)
  logpriorbeta1 <- dnorm(2,1,log=TRUE)
  logpriornu    <- dunif(phi,0.80,1.20,log=TRUE) 
  logprior <- logpriorbeta0+logpriorbeta1+
              logpriorphi+logpriornu
  like <- -(loglike/2) +logprior
  return(like)
}

##################
# Main M-H  loop #
##################

th <- c(0.01,0.01,0.01,0.01)

proposalfunction <- function(param,i,th){
  if (is.null(dim(param)[1])){
    sd <- th
    mu <- param
  }else{
    sd <- apply(param,2,sd)
    if (sum(sd) < 0.0001){
      sd <- th
    }else{
      sd <- sd
    }
    mu <- param[i,]
  }
  #alpha <- c(mu[1]^2/(th[1]*sd[1]))
  #beta  <- c(th[1]*sd[1]/mu[1])
  Yn <- c(#rgamma(1,shape=alpha[1],scale=beta[1]),
          runif(1,0.8,1),
          rnorm(1,mu[2],sd[2]),
          rnorm(1,mu[3],sd[3]),
          runif(1,0.8,1.20))
  return(Yn)
  #return(list(Yn,alpha,beta))
}

#metrop(f, startvalue, 10000)

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,npar))
  chain[1,] = startvalue
  for (i in 1:iterations){
    # iterations <- 10000;i<-1
    ## Decision
    proposal <- proposalfunction(chain[c(1:i),],i,th)
    #proposal <- proposal_all[[1]]
    #alpha <- proposal_all[[2]]
    #beta <- proposal_all[[3]]
    probab <- min(0,
    f(proposal) -
    #+ dgamma(chain[i,1],alpha[1],beta[1], log=TRUE)-
      f(chain[i,]) )
    #- dgamma(proposal[1],alpha[1],beta[1], log=TRUE))
    alphan <- exp(probab)
    if (log(runif(1)) <= probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
    if(i%%100==0){
    print(round(c(i, alphan, chain[i+1,]),4))
    }
  }
  return(chain)
}

print(datasetfile)

start_time <- Sys.time()

set.seed(100)
chain = run_metropolis_MCMC(startvalue, 2000)
burnIn = 500
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]));acceptance

end_time <- Sys.time()

print(end_time-start_time)

### Summary: #######################

png(filename=paste0("sim_res/plot",analysis,model,type,i,".png"))
par(mfrow = c(2,npar))
hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of phi", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),1]), col="green")
abline(v = 0.9, col="red" )
hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of beta0", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]), col="green")
abline(v = 0, col="red" )
hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of beta1", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]), col="green")
abline(v = 2, col="red" )
hist(chain[-(1:burnIn),4],nclass=30, main="Posterior of nu", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]), col="green")
abline(v = 1, col="red" )

plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of phi", )
abline(h = 0.9, col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of beta0", )
abline(h = 0, col="red" )
plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of beta1", )
abline(h = 2, col="red" )
plot(chain[-(1:burnIn),4], type = "l", xlab="True value = red line" , main = "Chain values of nu", )
abline(h = 1, col="red" )

dev.off()

save(chain, file=paste0("sim_res/chain",analysis,model,type,i,".Rdata"))
