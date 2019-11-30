#R code for M–H estimation 

##################
# L functions    #
##################

## log (Likelihood times prior)
f <- function(kappa,beta,taub,taue=1,a=1,b=0,m=2) {
loglike <- likelihoodKatzfuss(beta,kappa,sigma2,taue=1) #cambiar sigma2 por taub
logpriorbeta <- (0.5*log(taub)-0.5*taub(beta - m)^2) #1 beta
## incluir previas para taue y taub (según Demirhan et al)
logpriortaus <- dgamma(taub,shape=0.001, rate=1, lop=TRUE)  
# fijar k y theta para una previa de precisión Gamma.
logpriorkappa <- log(a - b) # fijar hiperparámetros de acuerdo al 1/rango
logprior <- logpriorbeta+logpriortaus+logpriorkappa
f <- loglike+logprior
return(f)
}

# para el miércoles
# plantear todo en términos de la precisión y no variancia. 
# Arreglar previa conjunta de las dos precisiones. 
# parametro rango y de precisión independencia en matern ** buscar referencias
# referencias apuntan a estudios de sensibilidad para ver si la independencia tiene sentido.

##################
#     Data       #
##################

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


for (t in 2:T) { 

## Candidate distribution for kappa,sigma2,beta
kappa.n <- kappa[t-1]+0.01*rnorm(1,0,1) ## proposals
sigma2.n <- sigma2[t-1]+0.01*rnorm(1,0,1) ## proposals
beta.n <- beta[t-1]+0.01*rnorm(1,0,1) ## proposals
tn <- f(kappa.n,sigma2.n,beta.n)
tc <- f(kappa[t-1],sigma2[t-1],beta[t-1])

## Decisión
if (log(u1[t]) <= tn-tc) 
  kappa[t] <- kappa.n; 
  prec[t] <- prec.n; 
  beta[t] <- beta.n  else
  {kappa[t] <- kappa[t-1]; 
    prec[t] <- prec[t-1]; 
    beta[t] <- beta[t-1]; 
    k1 <- k1+1 }

## Save values (tunning and/or sampling)

}

##################
#Summary & Stats #
##################

# posterior quantiles
quantile(kappa[1000:T], probs=c(.025,0.5,0.975))
quantile(prec[1000:T], probs=c(.025,0.5,0.975))
quantile(beta[1000:T], probs=c(.025,0.5,0.975))
# acceptance rates
1-k1/T