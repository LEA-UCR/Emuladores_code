install.packages("RcppArmadillo")
install.packages("Matrix")
install.packages("SparseM")
library(Matrix)
library(RcppArmadillo)
library(SparseM)
library(rbenchmark)

### Sparse Matrix example

data <- rnorm(1e6)
zero_index <- sample(1e6)[1:5e5]
data[zero_index] <- 0
mat <- matrix(data, ncol=1000)
mat[1:5,1:5]
mat_sparse <- Matrix(mat, sparse=TRUE)
class(mat_sparse)
dim(mat_sparse)

start_time <- Sys.time()
a<-solve(mat_sparse)
end_time <- Sys.time()
timeA1<-end_time-start_time

start_time <- Sys.time()
b<-mat_sparse%*%mat_sparse
end_time <- Sys.time()
timeA2<-end_time-start_time

c(timeA1,timeA2)

### Armadillo example


