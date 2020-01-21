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
  logpriortaue <- dunif(taue,1,30,log=TRUE) 
  #logpriortaue <- log(dinvgamma(taue,shape=5, scale=5))
  logpriortaub <- log(dinvgamma(taub,shape=5, scale=5))
  #logpriorphi <- dunif(phi,0.1,3,log=TRUE) 
  logprior <- logpriortaue+logpriortaub#+logpriorphi
  like <- -(loglike/2) +logprior
  return(like)
}

##################
# Main M-H  loop #
##################

Sn <- 0.001*diag(npar)
alphax <- 0.234

proposalfunction <- function(param, Un){
  Xn <- param
  Yn <- Xn+Sn%*%Un
  return(Yn)
}

#metrop(f, startvalue, nbatch = 1e2)

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,npar))
  chain[1,] = startvalue
  for (i in 1:iterations){
    # iterations <- 10000;i<-1
    ## Decision
    Un <- rnorm(npar,mean = chain[i,], sd= c(0.1,0.1))
    proposal = proposalfunction(chain[i,],Un)
    probab <- min(0,f(proposal) - f(chain[i,]))
    alphan <- exp(probab)
    if (log(runif(1)) <= probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
    ##Update according to Vihola 2012
    SSn <- Sn %*% (diag(npar)+(1/i)*c(alphan-alphax)*Un%*%t(Un)/
                     (norm(Un,'2')^2)) %*% t(Sn) 
    Sn <- t(chol(SSn))
    print(c(round(i,0), round(probab,4)))
  }
  return(chain)
}
chain = run_metropolis_MCMC(startvalue, 1000)
burnIn = 5
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
