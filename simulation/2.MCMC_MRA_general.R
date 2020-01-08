library(tictoc)
source("1.MRA_resolution_general.R")
source('../covariances.R')
source('likelihoodK_general.R')
library(MCMCpack)
library(truncdist)
library(invgamma)

#### Metropolis Hastings

# Variables:
js <- jb <- k1 <- 0
BurnIn <- 1000
TotIter <- 10000
u1 <- runif(TotIter)
AuxBurnIn <- 1
SaveResults <- list()


# initial values
nu <- 1.5
range <- 0.5
sigma2 <- 0.5
taue <- 20
betas <- 1
type='Exponential'


# storage space
nu.v <- rep(nu,TotIter)
range.v <- rep(range,TotIter)
taub.v <- rep(1/sigma2,TotIter)
taue.v <- rep(taue,TotIter)
beta.v <- rep(betas,TotIter)
tc.v <- rep(0,TotIter)
sd.v <- rep(0.01,TotIter)
##################
# L functions    #
##################

## log (-2*Likelihood times prior)
f <- function(nu,range,taub,taue=20,betas,type,a=0,b=1,m=2) {
  loglike <- likelihoodKatzfuss(nu,range,1/taub,taue,betas,type)
  #loglike <- likelihoodGaussian(nu,range,1/taub,taue,betas,type) #cambiar sigma2 por taub
  logpriorbeta <- -log(taub)+taub*(betas - m)^2 #1 beta
  ## incluir previas para taue y taub (según Demirhan et al). por ahora está en términos de la variancia, y no de la precisión
  logpriortaus <- -2*log(dinvgamma(taub,shape=0.001, scale=1))
  # fijar k y theta para una previa de precisión Gamma.
  logpriorrange <- -2*log(b - a) # fijar hiperparámetros de acuerdo al 1/rango
  logprior <- logpriorbeta+logpriortaus+logpriorrange
  like <- loglike+logprior
  return(like)
}

# pendiente para hacer
# plantear todo en términos de la precisión y no variancia. 
# Arreglar previa conjunta de las dos precisiones. 
# parametro rango y de precisión independencia:
# referencias apuntan a estudios de sensibilidad para ver si la independencia tiene sentido.

##################
# Main M-H  loop #
##################


#We need a starting value \theta^{(1)} and proposal density q(.,\theta_{i-1})
#For i in 2, .., B
#1. Draw a candidate \theta_{cand} ~q(.,\theta_{i-1})
#2. Compute r = \frac{p(y|\theta_{cand})p(\theta_{cand})q(\theta_{(i−1)};\theta_{cand}}
#                    {p(y|\theta_{(i−1)})p(\theta_{(i−1)})q(\theta_{cand};\theta_{(i−1)})}
#3. With probability set min(1,r), set \theta_i = \theta_{cand}.
#Otherwise, set \theta_i = \theta_{i-1}

tc.v[1] <- f(nu,range.v[1],taub.v[1],taue,beta.v[1],type)
for (t in 2:TotIter) { 
  
  ## Candidate distribution for kappa,sigma2,beta
  #range.n <- range.v[t-1]+0.01*rnorm(1,0,1) ## proposals
  range.n <- range.v[t-1]
  #range.n <- rtrunc(1000,spec = 'invgamma',a = 0,b = 1,shape=1,scale=1)  
  #taub.n <- taub.v[t-1]+0.01*rnorm(1,0,1) ## proposals
  taub.n <- taub.v[t-1]
  #if(is.na(sd(beta.v[1:(t-1)]))){
  #  sd.v[t] <- 0.01
  #}else{
  #  sd.v[t] <- sd(beta.v[1:(t-1)])
  #}
         
  beta.n <- beta.v[t-1]+sd.v[t]*rnorm(1,0,1) ## proposals
  tn <- f(nu,range.n,taub.n,taue,beta.n,type)
  
  ## Decisión
  decision <- min(0,tn-tc.v[t-1])
  show(decision)
  if (log(u1[t]) <= decision){
    range.v[t] <- range.n; 
    taub.v[t] <- taub.n; 
    beta.v[t] <- beta.n;
    tc.v[t] <- tn;
    k1 <- k1+1 }else
    {range.v[t] <- range.v[t-1]; 
    taub.v[t] <- taub.v[t-1]; 
    beta.v[t] <- beta.v[t-1];
    tc.v[t] <- tc.v[t-1]
    }
  
  ## Save values (tunning and/or sampling)
  print(k1)
}

##################
#Summary & Stats #
##################

#http://wlm.userweb.mwn.de/R/wlmRcoda.htm

# posterior quantiles
quantile(kappa.v[BurnIn:TotIter], probs=c(.025,0.5,0.975))
quantile(prec.v[BurnIn:TotIter], probs=c(.025,0.5,0.975))
quantile(betas.v[BurnIn:TotIter], probs=c(.025,0.5,0.975))
# acceptance rates
1-k1/TotIter
