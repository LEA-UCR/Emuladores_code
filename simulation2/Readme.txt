Plan de simulación

Generación de datos:

############################################
############################################

0.generate_data con cuatro modelos:
beta1s:: 2 options: Exponential and Mattern
GPerror:: 2 options: Exponential and Mattern
# Nugget effect generation
error <- rnorm(n*k, 0, sqrt(1/taue)) 
X <- rgamma(n*k,2,1) # X1
# Model beta(s)
y_BS <- beta0+(beta1)*scale(X) + (beta1s)*X + error 
X <- rep(1,n*k)
# Model GP # X2
y_GP <- beta0 + GPerror + error

OP1: beta0=0, beta1=2, taue = 10, range=0.9, 
sigma2=1, type="Exponential", X2
OP2: beta0=0, beta1=2, taue = 10, range=0.9, 
sigma2=1, type="Exponential", X1
OP3: beta0=0, beta1=2, taue = 10, range=0.9, 
sigma2=1, nu=1.5, type="Mattern", X2
OP4: beta0=0, beta1=2, taue = 10, range=0.9, 
sigma2=1, nu=1.5, type="Mattern", X1

############################################
############################################

## Tratamientos:
spBayes (MCMC_spBayes)
Likelihood (MCMC_LH)
FSA (MRA 1) (MCMC_MRA1)

############################################
############################################

## Response variable

MAE (for response Y)
MSE (for response Y)
Time (in seconds)



