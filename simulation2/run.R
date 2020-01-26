# if you are working local, setwd in simulation2 first!
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  i<-1
  type<-"Exponential"
  model<-"SVC"
  analysis<-"M1"
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

taub <- 1
taue <- 10

startvalue <- c(taub,taue)
N <- dim(hh)[1]
npar <- length(startvalue)

# fixed
phi <- 0.9
nu <- 1.5
beta0 <- 0
beta1 <- 2

##################
# L functions    #
##################

f <- function(param) {
  taub <- param[1]
  taue <- param[2]
  sigma2 <- 1/taub
  if(analysis=="M1"){
    loglike <- likelihoodGaussian(nu,phi,beta0,beta1,sigma2,taue,model,type)}
  if(analysis=="M2"){
    loglike <- likelihoodBanerjee(nu,phi,beta0,beta1,sigma2,taue,model,type)}
  else{loglike <- likelihoodFSA_Block(nu,phi,beta0,beta1,sigma2,taue,model,type)}
  #loglike <- likelihood(nu,phi,beta0,beta1,1/taub,taue,model,type)
  ## incluir previas para taue y taub (según Demirhan et al)

  logpriortaue <- (dgamma(taub,shape=0.5, scale=2, log=T))
  #logpriortaue <- log(dinvgamma(taue,shape=5, scale=5))
  logpriortaub <- (dgamma(taub,shape=5, scale=2, log=T))

  #logpriorphi <- dunif(phi,0.1,3,log=TRUE) 
  logprior <- logpriortaue+logpriortaub#+logpriorphi
  like <- -(loglike/2) +logprior
  return(like)
}

##################
# Main M-H  loop #
##################


th <- c(0.01,0.1)
alphax <- 0.234

proposalfunction <- function(param,i,th){
  if (is.null(dim(param)[1])){
    sd <- c(0.05,0.05)
    mu <- param
  }else{
    sd <- apply(param,2,sd)
    mu <- param[i,]
  }
  alpha <- c(mu[1]^2/(th[1]*sd[1]),mu[2]^2/(th[2]*sd[2]))
  beta  <- c(th[1]*sd[1]/mu[1],th[2]*sd[2]/mu[2])
  Yn <- c(rgamma(1,shape=alpha[1],scale=beta[1]),
          rgamma(1,shape=alpha[2],scale=beta[2]))
  return(list(Yn,alpha,beta))
}

#metrop(f, startvalue, nbatch = 1e2)

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,npar))
  chain[1,] = startvalue
  for (i in 1:iterations){
    # iterations <- 10000;i<-1
    ## Decision
    proposal_all <- proposalfunction(chain[c(1:i),],i,th)
    proposal <- proposal_all[[1]]
    alpha <- proposal_all[[2]]
    beta <- proposal_all[[3]]
    probab <- min(0,
    f(proposal) + dgamma(chain[i,1],alpha[1],beta[1], log=TRUE)
        + dgamma(chain[i,2],alpha[2],beta[2], log=TRUE)
  -f(chain[i,]) - dgamma(proposal[1],alpha[1],beta[1], log=TRUE)
        - dgamma(proposal[2],alpha[2],beta[2], log=TRUE))
    alphan <- exp(probab)
    if (log(runif(1)) <= probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }

    print(c(round(i,0), round(alphan,4), round(chain[i+1,],4)))
  }
  return(chain)
}
chain = run_metropolis_MCMC(startvalue, 40000)
burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]));acceptance

### Summary: #######################

png(filename=paste0("sim_res/plot",analysis,model,type,i,".png"))
par(mfrow = c(2,npar))
# hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of phi", xlab="True value = red line")
# abline(v = mean(chain[-(1:burnIn),2]))
# abline(v = 0.9, col="red" )
hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of taub", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]))
abline(v = 1, col="red" )
hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of taue", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]) )
abline(v = 10, col="red" )
# plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of phi", )
# abline(h = 0.9, col="red" )
plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of taub", )
abline(h = 1, col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of taue", )
abline(h = 10, col="red" )
dev.off()

save(chain, file=paste0("sim_res/chain",analysis,model,type,i,".Rdata"))
