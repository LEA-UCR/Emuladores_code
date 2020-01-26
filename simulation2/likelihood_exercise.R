# if you are working local, setwd in simulation2 first!
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  i<-1
  type<-"Exponential"
  model<-"SVC"
  analysis<-"M1"
  datasetfile=paste0("sim_data/dataset",model,type,i,".Rdata")
} else {
  i<-args[1]
  type<-args[2]
  model<-args[3]
  analysis<-args[4]
  datasetfile=paste0("sim_data/dataset",model,type,i,".Rdata")
}

# i<-1:100
# type<-'Exponential', "Matern"
# model<-'SVC', "SVI"
# analysis<-"M1: likelihood", "M2: Banerjee", "M3: FSA"

source("1.MRA_resolution_general.R")
source('covariances.R')
source('likelihoodK_general.R')
library(MCMCpack)
library(truncdist)
library(invgamma)

aa<-gen_resolution(datasetfile)
bordes<-aa[[1]];indicesW<-aa[[2]];knotsMRA<-aa[[3]]
nn<-aa[[4]];hh<-aa[[5]]

#### Likelihood tests

# initial values

# taub <- 1
# taue <- 10

taub<-seq(0.5,4,0.25)
taue<-seq(1,30,1)

parameters<-expand.grid(taub,taue)

# fixed
phi <- 0.9
nu <- 1.5
beta0 <- 0
beta1 <- 2

################################
#### evaluate L functions   ####
################################

M1<-function(param){
  taue <- param[,1]
  taub <- param[,2]
  sigma2 <- 1/taub
  m2logv <-likelihoodGaussian(nu,phi,beta0,beta1,sigma2,taue,model,type)
  return(m2logv)}
M2 <- function(param){
  taue <- param[,1]
  taub <- param[,2]
  sigma2 <- 1/taub
  m2logv<-likelihoodBanerjee(nu,phi,beta0,beta1,sigma2,taue,model,type)
  return(m2logv)}
M3 <- function(param){
  taue <- param[,1]
  taub <- param[,2]
  sigma2 <- 1/taub
  m2logv<-likelihoodFSA_Block(nu,phi,beta0,beta1,sigma2,taue,model,type)
  return(m2logv)}


resM1<-unlist(lapply(1:dim(parameters)[1],
                     function(i)M1(parameters[i,])))
resM2<-unlist(lapply(1:dim(parameters)[1],
                     function(i)M2(parameters[i,])))
resM3<-unlist(lapply(1:dim(parameters)[1],
                     function(i)M3(parameters[i,])))


results<-as_tibble(cbind(taue=parameters[,1],taub=parameters[,2],
                         M1=-resM1/2,M2=-resM2/2,M3=-resM3/2))

par(mfrow=c(1,3))
aa<-matrix(results$M1,length(unique(results$taue)),
           length(unique(results$taub)))
image.plot(x = unique(results$taub), y = unique(results$taue),
             z = t(aa), main="Maximizing Gaussian Likelihood",
           ylab=expression(tau^2*epsilon),
           xlab=expression(tau^2*beta))
abline (v=10,h=1)

aa<-matrix(results$M2,length(unique(results$taue)),
           length(unique(results$taub)))
image.plot(x = unique(results$taub), y = unique(results$taue),

           z = t(aa), main="Maximizing Banerjee Likelihood",
           ylab=expression(tau^2*epsilon),
           xlab=expression(tau^2*beta))

abline (v=10,h=1)

aa<-matrix(results$M3,length(unique(results$taue)),
           length(unique(results$taub)))
aa<-ifelse(aa==Inf,0,aa)
image.plot(x = unique(results$taub), y = unique(results$taue),

           z = t(aa), main="Maximizing FSA Likelihood",
           ylab=expression(tau^2*epsilon),
           xlab=expression(tau^2*beta))
abline (v=10,h=1)

## time
library(tictoc)
tic()
M1(parameters[5,])
toc()
tic()
M2(parameters[5,])
toc()
tic()
M3(parameters[5,])
toc()

