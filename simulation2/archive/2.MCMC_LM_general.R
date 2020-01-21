library(tictoc)
source("simulation2/1.MRA_resolution_general.R")
source('covariances.R')
source('simulation2/likelihoodK_general.R')
library(MCMCpack)
library(truncdist)
library(invgamma)

### ############################ ###
###      Metropolis Hastings     ###
### ############################ ###

##################
# L functions    #
##################

likelihood <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  x <- hh$X
  y <- hh$Y_LM
  pred = a*scale(x) + b
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)   
}

prior <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  aprior = dunif(a, min=0, max=10, log = T)
  bprior = dnorm(b, sd = 5, log = T)
  sdprior = dunif(sd, min=0, max=30, log = T)
  return(aprior+bprior+sdprior)
}

posterior <- function(param){
  return (likelihood(param) + prior(param))
}

##################
# Main M-H  loop #
##################

proposalfunction <- function(param){
  return(rnorm(3,mean = param, sd= c(0.1,0.5,0.3)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

startvalue <- c(2,0,10)
chain = run_metropolis_MCMC(startvalue, 100000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]));acceptance

##################
#Summary & Stats #
##################

### Summary: #######################

par(mfrow = c(2,3))
hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of beta1", xlab="True value = red line")
abline(v = mean(chain[-(1:BurnIn),2]))
abline(v = 2, col="red" )
hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of beta0", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]))
abline(v = 0, col="red" )
hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of sigma_e", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]) )
abline(v = sqrt(10), col="red" )
plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of beta1", )
abline(h = 2, col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of beta0", )
abline(h = 0, col="red" )
plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sigma_e", )
abline(h = sqrt(10), col="red" )

summary(lm(hh$Y_LM~hh$X_scaled))$sigma
summary(lm(hh$Y_LM~hh$X_scaled))$coef
