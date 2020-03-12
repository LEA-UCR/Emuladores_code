rm(list=ls())
library(readr)
library(coda)

real_values <- c(0.9,0,2,1)
i<-1:20
type<-c("Matern","Exponential")
model<-c("SVI","SVC")
analysis<-c("M1","M2")

extract_all <- function(analysis, model, type, i){
  file1 <- paste0("sim_res/chain",analysis,model,
                  type,i,".Rdata")
  file2 <- paste0("sim_res/output_",i,"_",type,
                  "_",model,"_",analysis,".txt")
  load(file1)
  
  vars <- which(abs(geweke.diag(chain, frac1=0.1, frac2=0.5)$z)<
                  qnorm(0.975))
  
  median_c <- apply(chain,2,median)
  mean_c <- apply(chain,2,mean)
  sd_c <- apply(chain,2,sd)
  n_c <- apply(chain,2,length)-1
  descriptives <- cbind(median_c,mean_c,sd_c,n_c)
  
  bias<-apply(chain,2,mean)-real_values
  mse_func <- function(x){sum((x-mean(x))^2/length(x))}
  mse <- apply(chain,2,mse_func)
  mad_func <- function(x){median(abs(x-median(x)))}
  mad <- apply(chain,2,mad_func)
  time <- parse_number(as.character(read.delim(file2)[303,]))
  mcmc_desc <- cbind(bias,mse,mad)
  
  return(list(descriptives,mcmc_desc,time, converged=vars))
}


analysis<-c("M1","M2","M3")
model<-c("SVI","SVC")
type<-c("Matern","Exponential")
i<-1:20

extract_all(analysis[1],model[2],type[2],i[1])[[3]]
extract_all(analysis[1],model[2],type[2],i[2])[[3]]
extract_all(analysis[1],model[2],type[2],i[3])[[3]]
extract_all(analysis[1],model[2],type[2],i[4])[[3]]

extract_all(analysis[1],model[2],type[1],i[1])[[3]]

extract_all(analysis[1],model[1],type[2],i[1])[[3]]
extract_all(analysis[1],model[1],type[2],i[4])[[3]]

extract_all(analysis[1],model[1],type[1],i[1])[[3]]

extract_all(analysis[2],model[2],type[2],i[1])[[3]]
extract_all(analysis[2],model[2],type[2],i[2])[[3]]
extract_all(analysis[2],model[2],type[2],i[3])[[3]]
extract_all(analysis[2],model[2],type[2],i[4])[[3]]

extract_all(analysis[2],model[2],type[1],i[10])[[3]]

#extract_all(analysis[3],model[2],type[1],i[1])[[3]]
#extract_all(analysis[3],model[1],type[2],i[1])[[3]]
#extract_all(analysis[3],model[1],type[1],i[1])[[3]]