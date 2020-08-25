# if you are working local, setwd in the working folder first!
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  variable_narccap <- 'Temp'
  type<-"Exponential"
  model<-"SVC"
  analysis<-"M3"
  datasetfile=paste0("data_narccap/dataset",
                     variable_narccap,".Rdata")
} else {
  variable_narccap <-args[1]
  type<-args[2]
  model<-args[3]
  analysis<-args[4]
  datasetfile=paste0("data_narccap/dataset",
                     variable_narccap,".Rdata")
}

# i<-1:100
# type<-'Exponential', "Matern"
# model<-'SVC', "SVI"
# analysis<-"M1: likelihood", "M2: FSA", "M3: MRA2"

source("1.MRA_resolution_general.R")
source('covariances.R')
source('likelihoodK_general.R')
library(MCMCpack)
library(truncdist)
library(invgamma)
library(Matrix)
library(mvtnorm)

nCov_f <- 6
nCov_v <- 5

aa<-gen_resolution(datasetfile,nCov_v)
Tm <- aa[[1]]
TSm <- aa[[2]]
Qlist <- aa[[3]]
nn<-aa[[4]]
hh <- Qlist[[nn+2]]

#### Metropolis Hastings

# initial values

phi <- 0.9
nu <- 1
beta <- c(0,rep(1,nCov_f-1))

startvalue <- c(phi,nu,beta)
N <- dim(hh)[1]
npar <- length(startvalue)

# fixed
taub <- 1
taue <- 1
sigma2 <- 1/taub

A <- sigma2*diag(nCov_v)
types <- rep('Exponential',nCov_v)

#Data
Y <- hh$Y
X <- data.matrix(st_drop_geometry(hh %>% mutate(interc = 1) %>% 
                      dplyr::select(interc,TREFHT,OMEGA,PSL,U,V)))
XR <- scale(X[,c(2,3,4,5,6)])

##################
# L functions    #
##################

f <- function(param) {
  phi   <- param[1]
  nu    <- param[2]
  beta <- param[-c(1,2)]
  if (analysis=="M1"){
    loglike <- likelihoodGaussian(nu,phi,beta0,
                  beta1,sigma2,taue,model,type)
    }else{
      if (analysis=="M2"){
    loglike <- likelihoodFSA_Block(nu,phi,beta0,
                  beta1,sigma2,taue,model,type)
      }else {
        MRA_num <- nn
        loglike <- likelihoodMRA(nu,phi,beta,A,nCov_v,
                                 taue,model,type, MRA_num,Y,X,XR)}}
  logpriorphi   <- dunif(phi,0.80,1.00,log=TRUE) 
  logpriornu    <- dunif(nu,0.80,1.20,log=TRUE) 
  logpriorbeta  <- dmvnorm(beta,mean = c(0,rep(1,length(beta)-1)),
                           sigma = 0.5*diag(length(beta)),log = T)
  logprior <- logpriorphi+logpriornu+logpriorbeta
  like <- -(loglike/2) +logprior
  return(as.numeric(like))
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
          runif(1,0.8,1.20),
          rmvnorm(1,mean = mu[-c(1,2)],sigma = diag(mu[-c(1,2)])))
  return(Yn)
  #return(list(Yn,alpha,beta))
}

#metrop(f, startvalue, 10000)

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,npar))
  chain[1,] = startvalue
  fchain <- f(chain[1,])
  for (i in 1:iterations){
    show(i)
    # iterations <- 10000;i<-1
    ## Decision
    proposal <- proposalfunction(chain[c(1:i),],i,th)
    #proposal <- proposal_all[[1]]
    #alpha <- proposal_all[[2]]
    #beta <- proposal_all[[3]]
    fprop <- f(proposal)
    probab <- min(0,
    fprop -
    #+ dgamma(chain[i,1],alpha[1],beta[1], log=TRUE)-
      fchain)
    #- dgamma(proposal[1],alpha[1],beta[1], log=TRUE))
    alphan <- exp(probab)
    if (log(runif(1)) <= probab){
      chain[i+1,] = proposal
      fchain <- fprop
    }else{
      chain[i+1,] = chain[i,]
    }
    if(i%%100==0){
    print(round(c(i, alphan, chain[i+1,]),4))
    }
  }
  return(chain)
}

print(paste("Model =",analysis,"/ Data =",datasetfile))

start_time <- Sys.time()

set.seed(100)
chain = run_metropolis_MCMC(startvalue, 1000)
burnIn = 200
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]));acceptance

end_time <- Sys.time()

print(end_time-start_time)

### Summary: #######################

png(filename=paste0("sim_res/plot",analysis,model,type,'NARCCAP',".png"))
par(mfrow = c(2,npar))
hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of phi", 
     xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),1]), col="green")
abline(v = 0.9, col="red" )
hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of beta0", 
     xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]), col="green")
abline(v = 0, col="red" )
hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of beta1", 
     xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]), col="green")
abline(v = 2, col="red" )
hist(chain[-(1:burnIn),4],nclass=30, main="Posterior of nu", 
     xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),4]), col="green")
abline(v = 1, col="red" )

plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , 
     main = "Chain values of phi", )
abline(h = 0.9, col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , 
     main = "Chain values of beta0", )
abline(h = 0, col="red" )
plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , 
     main = "Chain values of beta1", )
abline(h = 2, col="red" )
plot(chain[-(1:burnIn),4], type = "l", xlab="True value = red line" , 
     main = "Chain values of nu", )
abline(h = 1, col="red" )

dev.off()

save(chain, file=paste0("sim_res/chain",
    analysis,model,type,'NARCCAP',".Rdata"))
