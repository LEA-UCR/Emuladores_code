set.seed(1000)

# Genero datos
N <- 1000
X <- matrix(rnorm(N*N),N,N)
y <- X%*%rep(8,N)+rnorm(N,0,1)
MCMCBetasI <- rep(1,N)

# Hiperparámetros (fijos)
sigma2_ <- 1
var <- 1
V <- diag(N)*var
b <- 1
a <- 0
vi <- 1
mi <- 0

# otras cantidades necesarias:
js <- jb <- 0
BurnIn <- 1000
TotIter <- 10000
AuxBurnIn <- 1
SaveResults <- list()

#https://brunaw.com/phd/mcmc/report.pdf

for(i in 1:TotIter){
  ## Proposal distribution for sigma2

  sigma2c = rgamma(1, (sigma2_^2)/var, scale=1/(sigma2_ / var))
  # Computing the ratio in the Metropolis-Hastings
  ratio_sigma2 = (((-N/2) * log(sigma2c) - 0.5/sigma2c *
   (t(y - (X%*%MCMCBetasI))%*%(y-(X%*%MCMCBetasI))) - 
     log(b - a) + 1/sigma2_ + 
     dgamma(sigma2_, (sigma2c^2)/var, scale = 1/(sigma2c / var))) -
     ((-N/2) * log(sigma2_) - 0.5/sigma2_ * 
        (t(y - (X%*%MCMCBetasI))%*%(y-(X%*%MCMCBetasI))) -
        log(b - a) + 1/sigma2c + 
        dgamma(sigma2c, (sigma2_^2)/var, scale = 1/(sigma2_ / var))))
  # Accept/Reject step for sigma2
  if(runif(1) < min(1, exp(ratio_sigma2)))
  {
    sigma2_ = sigma2c
    js = js +1
  }
  
  ## Proposal distribution for beta
  
  MCMCBetasC <- sapply(MCMCBetasI, function(x)rnorm(1,x,1))
  
  # Computing the ratio:
  ratio_beta <- ((-0.5/sigma2_* (t(y - (X %*% MCMCBetasC)) %*%
    (y - (X %*% MCMCBetasC)) - 0.5/vi*sum(abs(MCMCBetasC - mi))^2)) -
      (-0.5/sigma2_* (t(y - (X %*% MCMCBetasI)) %*%
  (y - (X %*% MCMCBetasI)) - 0.5/vi*sum(abs(MCMCBetasI - mi))^2)))
  # Accept/Reject step for beta
  if(runif(1) < min(1, exp(ratio_beta)))
  {
    MCMCBetasI <- MCMCBetasC
    jb <- jb+1
  }
  
  ## Leave burned runs out
  if (i > BurnIn){
    # Saving results
    SaveResults[[AuxBurnIn]] <- c(AuxBurnIn, MCMCBetasI, sigma2_)
    AuxBurnIn <- AuxBurnIn + 1
  }
  print(i)
}

### See results: use coda package
### 


library(coda)
# http://www.math.kit.edu/stoch/lehre/abib2010w/media/coda.pdf
samples_beta <- as.mcmc(t(matrix(unlist(SaveResults),nrow=1002))[,2:1001])
samples_sigma <- as.mcmc(t(matrix(unlist(SaveResults),nrow=1002))[,1002])
plot(samples_sigma)
matplot(t(samples_beta), type="l")
raftery.diag(samples)
effectiveSize(samples)
rejectionRate(samples)
