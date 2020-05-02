rm(list=ls())
library(readr)
library(coda)
library(purrr)
library(tidyverse)

#setwd("~/Dropbox/Emuladores/Emuladores_code/simulation3")

real_values <- c(0.9,0,2,1)
i<-1:10
type<-c("Matern","Exponential")
model<-c("SVI","SVC")
analysis<-c("M1","M2","M3")

extract_all <- function(analy, mod, typ, ind){
  file1 <- paste0("sim_res/chain",analy,mod,
                  typ,ind,".Rdata")
  #file2 <- paste0("sim_res/output_",ind,"_",typ,
  #                "_",mod,"_",analy,".txt")
  load(file1)
  
  vars <- abs(geweke.diag(chain, frac1=0.1, frac2=0.5)$z)<
                  qnorm(0.995)
  
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
 # if(analy=="M3"){
 # time <- parse_number(as.character(read.delim(file2)[[1]][100103]))}
 #   else{
 # time <- parse_number(as.character(read.delim(file2)[303,]))}
  mcmc_desc <- cbind(bias,mse,mad)
  
  return(descript=cbind(descriptives,mcmc_desc,converged=vars))
}


analysis<-c("M1","M2","M3")
model<-c("SVI","SVC")
type<-c("Matern","Exponential")
i<-1:10

OP1<-map(1:10,function(j){extract_all(analysis[1],model[1],type[1],j)})
OP2<-map(1:10,function(j){extract_all(analysis[1],model[1],type[2],j)})
OP3<-map(1:10,function(j){extract_all(analysis[1],model[2],type[1],j)})
OP4<-map(1:10,function(j){extract_all(analysis[1],model[2],type[2],j)})
OP5<-map(1:10,function(j){extract_all(analysis[2],model[1],type[1],j)})
OP6<-map(1:10,function(j){extract_all(analysis[2],model[1],type[2],j)})
OP7<-map(1:10,function(j){extract_all(analysis[2],model[2],type[1],j)})
OP8<-map(1:10,function(j){extract_all(analysis[2],model[2],type[2],j)})
OP9<-map(1:10,function(j){extract_all(analysis[3],model[1],type[1],j)})
OP10<-map(1:10,function(j){extract_all(analysis[3],model[1],type[2],j)})
OP11<-map(1:10,function(j){extract_all(analysis[3],model[2],type[1],j)})
OP12<-map(1:10,function(j){extract_all(analysis[3],model[2],type[2],j)})


datos<-as_tibble(rbind(cbind(do.call(rbind,OP1),setting=rep(1,40)),
      cbind(do.call(rbind,OP2),setting=rep(2,40)),
      cbind(do.call(rbind,OP3),setting=rep(3,40)),
      cbind(do.call(rbind,OP4),setting=rep(4,40)),
      cbind(do.call(rbind,OP5),setting=rep(5,40)),
      cbind(do.call(rbind,OP6),setting=rep(6,40)),
      cbind(do.call(rbind,OP7),setting=rep(7,40)),
      cbind(do.call(rbind,OP8),setting=rep(8,40)),
      cbind(do.call(rbind,OP9),setting=rep(9,40)),
      cbind(do.call(rbind,OP10),setting=rep(10,40)),
      cbind(do.call(rbind,OP11),setting=rep(11,40)),
      cbind(do.call(rbind,OP12),setting=rep(12,40))))

datos_f<-datos %>% mutate(var=rep(1:4,120), rep=rep(rep(1:10,each=4),12)) #%>% 
#  mutate(time_10x = time_chain*rep(c(rep(1/3,8),rep(1,4)),each=40))
summary(datos_f)
unique(datos_f %>% filter(converged==0) %>% select(var,rep,setting))

# time for 10000 iterations for each setting:
#datos_f %>% group_by(setting) %>% 
#  select(time_10x,setting) %>% 
#  summarise(mean_time=mean(time_10x), 
#            sd_time=sd(time_10x)) 
  
a<-datos_f %>% group_by(var,setting) %>% 
  select(bias,mse,mad,setting) %>% 
  summarise(mean_bias=mean(bias), 
            sd_bias=sd(bias),
            mean_mse=mean(mse), 
            sd_mse=sd(mse),
            mean_mad=mean(mad), 
            sd_mad=sd(mad)) 
a
print(a, n=48)
