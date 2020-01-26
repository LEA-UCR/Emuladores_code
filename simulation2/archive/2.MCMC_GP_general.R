datasetfile <- "simulation2/sim_data/datasetSVCExponential1.Rdata"
method <- 'simulation2/likelihoodK_general.R'

source("simulation2/1.MRA_resolution_general.R")
source('covariances.R')
source(method)
library(MCMCpack)
library(truncdist)
library(invgamma)

aa<-gen_resolution(datasetfile)
bordes<-aa[[1]];indicesW<-aa[[2]];knotsMRA<-aa[[3]]
nn<-aa[[4]];hh<-aa[[5]]

### ############################ ###
###      Metropolis Hastings     ###
### ############################ ###

# parameters: beta0=0, beta1=2.
# Spatial par: taue=10, range=0.9, 
# Spatial par: taub=sigma2=1,
# Fixed: nu=1 (fixed), type='Exponential'(fixed)

# c(beta1,beta0,taue,range,sigma2)
startvalue <- #c(2,0,
  c(10,0.9)#,1)
N <- dim(hh)[1]
npar <- length(startvalue)

##################
# L functions    #
##################

likelihood <- function(param){
  beta1 = 2#param[1]
  beta0 = 0#param[2]
  taue = param[1]
  range = param[2]
  sigma2 = 1#param[5]
  x <- hh$X_scaled
  y <- hh$Y_LM
  pred = beta1*x + beta0
  Sigma_e <- cExpMat(hh,hh,type='Exponential',range,(1/taue),nu=1)+sigma2*diag(N)
  Sigmainv <- chol2inv(Sigma_e)
  logv <- -200*log(2*pi)+(0.5)*log(det(Sigmainv))-(0.5)*(t(y-pred)%*%Sigma_e%*%(y-pred))
  return(logv)   
}

prior <- function(param){
  #beta1 = param[1]
  #beta0 = param[2]
  taue = param[1]
  range = param[2]
  #sigma2 = param[5]
  #aprior = dunif(beta1, min=1.5, max=2.5, log = T)
  #bprior = dnorm(beta0, sd = 0.01, log = T)
  rprior = dunif(range, min=0.5, max=1.1, log = T) #dgamma(range, 5, 5, log=TRUE)
  tprior = log(dinvgamma(taue, 5, 40))
  #s2prior = log(dinvgamma(sigma2, 5, 4))
  return(#aprior+bprior+
    rprior+tprior)#+s2prior)
}

posterior <- function(param){
  return (likelihood(param) + prior(param))
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

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,npar))
  chain[1,] = startvalue
  for (i in 1:iterations){
    # iterations <- 10000;i<-1
    ## Decisión
    Un <- rnorm(npar,mean = chain[i,], sd= c(#0.0001,0.00001,
      0.1,0.01))#,0.01))
    proposal = proposalfunction(chain[i,],Un)
    probab <- posterior(proposal) - posterior(chain[i,])
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

chain = run_metropolis_MCMC(startvalue, 10000)

burnIn = 5
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]));acceptance

##################
#Summary & Stats #
##################

### Summary: #######################

par(mfrow = c(2,npar))

# hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of beta1", xlab="True value = red line")
# abline(v = mean(chain[-(1:burnIn),2]))
# abline(v = 2, col="red" )
# hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of beta0", xlab="True value = red line")
# abline(v = mean(chain[-(1:burnIn),2]))
# abline(v = 0, col="red" )
hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of tau_e", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]) )
abline(v = 10, col="red" )
hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of range", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]) )
abline(v = 0.9, col="red" )
# hist(chain[-(1:burnIn),5],nclass=30, main="Posterior of sigma2", xlab="True value = red line")
# abline(v = mean(chain[-(1:burnIn),2]) )
# abline(v = 1, col="red" )
# plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of beta1", )
# abline(h = 2, col="red" )
# plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of beta0", )
# abline(h = 0, col="red" )
plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of tau_e", )
abline(h = 10, col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of range", )
abline(h = 0.9, col="red" )
# plot(chain[-(1:burnIn),5], type = "l", xlab="True value = red line" , main = "Chain values of sigma2", )
# abline(h = 1, col="red" )

# c(beta1,beta0,taue,range,sigma2)

### ############################ ###
###  With spBayes (same model)   ###
### ############################ ###

library(spBayes)

# number of MCMC iterations, used in spLM()
n.samples <- 10000
# starting values for parameters in spLM()
starting <- list("phi"=0.9, "sigma.sq"=1, "tau.sq"=10)
# variances for the Metropolis sampler in spLM()
tuning <- list("phi"=0.001, "sigma.sq"=0.01, "tau.sq"=0.1)
# priors for parameters in spLM(): betas are multivariate Normal
#priors.1 <- list("beta.Norm"=list(rep(0,2), diag(1000,2)),
#                 "phi.Unif"=c(0.01, 0.5), "sigma.sq.IG"=c(5, 5),
#                 "tau.sq.IG"=c(5, 5))
# priors for parameters in spLM(): betas are flat
priors.2 <- list("beta.Flat", "phi.Unif"=c(0.5, 2),
                 "sigma.sq.IG"=c(5, 5), "tau.sq.IG"=c(5, 50))
# function for spatial dependence structure in spLM()
cov.model <- "exponential"
# interval for seeing progress of the sampler in spLM()
n.report <- 500
coords <- expand.grid(seq(0,1,length.out = 20),
                      seq(0,1,length.out = 20))
# model with first set of priors
m.1 <- spLM(hh$Y_LM~hh$X_scaled-1, coords=as.matrix(coords), starting=starting,
            tuning=tuning, priors=priors.2, cov.model=cov.model,
            n.samples=n.samples, n.report=n.report)

par(mfrow=c(2,2))
ts.plot(m.1$p.theta.samples[,1],main="sigma sq",ylab="",
        xlim=c(100,nrow(m.1$p.theta.samples)),ylim=c(0,3))
ts.plot(m.1$p.theta.samples[,2],main="tau sq",ylab="",
        xlim=c(100,nrow(m.1$p.theta.samples)),ylim=c(0,3))
ts.plot(m.1$p.theta.samples[,3],main="phi",ylab="",
        xlim=c(50,nrow(m.1$p.theta.samples)))

