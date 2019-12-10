#R code for M–H estimation 

##################
#     Data       #
##################

library(MCMCpack)
source(file="0.generate_data.R")
source(file="2.MCMC_MRA.R")

N <- dim(regionalpoints)[1]
V <- diag(N)*sigma2[1]
b <- 1
a <- 0
vi <- 1
mi <- 0

# otras cantidades necesarias:
js <- jb <- 0
BurnIn <- 1000
TotIter <- 10000
u1 <- runif(TotIter)
AuxBurnIn <- 1
SaveResults <- list()

# initial values

kappa <- 1.5
prec <- 4
betas <- 1
sigma2 <- 1/4

kappa.v <- rep(kappa,TotIter)
prec.v <- rep(prec,TotIter)
betas.v <- rep(betas,TotIter)


##################
# L functions    #
##################

## log (Likelihood times prior)
f <- function(kappa,betas,taub,taue=1,a=1,b=0,m=2) {
  loglike <- likelihoodKatzfuss(betas,kappa,sigma2,taue=1) #cambiar sigma2 por taub
  logpriorbeta <- (0.5*log(taub)-0.5*taub*(betas - m)^2) #1 beta
  ## incluir previas para taue y taub (según Demirhan et al). por ahora está en términos de la variancia, y no de la precisión
  logpriortaus <- log(dinvgamma(taub,shape=0.001, scale=1))
  # fijar k y theta para una previa de precisión Gamma.
  logpriorkappa <- log(a - b) # fijar hiperparámetros de acuerdo al 1/rango
  logprior <- logpriorbeta+logpriortaus+logpriorkappa
  like <- loglike+logprior
  return(like)
}

# pendiente para hacer
# plantear todo en términos de la precisión y no variancia. 
# Arreglar previa conjunta de las dos precisiones. 
# parametro rango y de precisión independencia en matern ** buscar referencias
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

for (t in 2:TotIter) { 

## Candidate distribution for kappa,sigma2,beta
kappa.n <- kappa.v[t-1]+0.01*rnorm(1,0,1) ## proposals
prec.n <- prec.v[t-1]+0.01*rnorm(1,0,1) ## proposals
beta.n <- betas.v[t-1]+0.01*rnorm(1,0,1) ## proposals
tn <- f(kappa.n,beta.n,prec.n)
tc <- f(kappa[t-1],betas[t-1],prec[t-1])

## Decisión
if (log(u1[t]) <= tn-tc){
  kappa.v[t] <- kappa.n; 
  prec.v[t] <- prec.n; 
  betas.v[t] <- beta.n}else
  {kappa.v[t] <- kappa[t-1]; 
    prec.v[t] <- prec[t-1]; 
    beta.v[t] <- beta[t-1]; 
    k1 <- k1+1 }

## Save values (tunning and/or sampling)

}

##################
#Summary & Stats #
##################

#http://wlm.userweb.mwn.de/R/wlmRcoda.htm

# posterior quantiles
quantile(kappa.v[1000:T], probs=c(.025,0.5,0.975))
quantile(prec.v[1000:T], probs=c(.025,0.5,0.975))
quantile(beta.v[1000:T], probs=c(.025,0.5,0.975))
# acceptance rates
1-k1/T