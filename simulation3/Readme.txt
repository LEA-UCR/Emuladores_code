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
# Model GP 
X <- rep(1,n*k)
y_GP <- beta1*X + (beta1s)*X + error

OP1: beta0=0, beta1=2, nu = 1, range=0.9, 
sigma2 y taub fijos, type="Exponential", X2
OP2: beta0=0, beta1=2, nu = 1, range=0.9, 
sigma2 y taue fijos, type="Exponential", X1
OP3: beta0=0, beta1=2, nu = 1, range=0.9, 
sigma2 y taue fijos, type="Mattern", X2
OP4: beta0=0, beta1=2, nu = 1, range=0.9, 
sigma2 y taue fijos, type="Mattern", X1

############################################
############################################

## Tratamientos:

Likelihood (MCMC_LH)
FSA (MRA 1) (MCMC_MRA1)

############################################
############################################

## Response variable

Average Bias for beta0, beta1, range, and nu per each of the 8 options.
MAE (for beta0, beta1, range, and nu.) with respect to real value
MAD (for beta0, beta1, range, and nu.) calculate for each chain, then average.
Time (in seconds)



