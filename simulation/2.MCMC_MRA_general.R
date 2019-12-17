library(tictoc)
source("1.MRA_resolution_general.R")
source('../covariances.R')
source('likelihoodK_general.R')


# values used to generate the data
nu <- 1.5
range <- 4
sigma2 <- 1
taue <- 20
betas <- 1
type='Exponential'



tic()
res <- likelihoodKatzfuss(betas,kappa,sigma2,taue)
toc()

#dj_1,...,j_M    *X

#uj_1,...,j_M    *X

#Atilde k,l j_1,...,j_M     *X
#omegatilde k,l j_1,...,j_M      *X

#A k,l j_1,...,j_M-1
#omega k j_1,...,j_M-1

#Ktilde j_1,...,j_M-1

#dj_1,...,j_M-1
#uj_1,...,j_M-1

#Atilde k,l j_1,...,j_M-1
#omegatilde k j_1,...,j_M-1



