# install.packages("RcppArmadillo")
# install.packages("rbenchmark")
# install.packages("inline")
# install.packages("Matrix")

library(Matrix)
library(RcppArmadillo)
library(rbenchmark)
library(inline)

## References: https://cran.r-project.org/web/packages/RcppArmadillo/vignettes/RcppArmadillo-intro.pdf
## https://eddelbuettel.github.io/pinp/Rcpp-introduction.pdf
## https://cran.r-project.org/web/packages/inline/index.html

### Sparse Matrix example

data <- rnorm(1e6)
zero_index <- sample(1e6)[1:5e5]
data[zero_index] <- 0
mat <- matrix(data, ncol=1000)
mat[1:5,1:5]
mat_sparse <- Matrix(mat, sparse=TRUE)
class(mat_sparse)
dim(mat_sparse)

# start_time <- Sys.time()
# a<-solve(mat_sparse)
# end_time <- Sys.time()
# timeA1<-end_time-start_time
# 
# start_time <- Sys.time()
# b<-mat_sparse%*%mat_sparse
# end_time <- Sys.time()
# timeA2<-end_time-start_time
# 
# c(timeA1,timeA2)

### Armadillo example

mat[1:5,1:5]
g1 <- cxxfunction(signature (vs = "numeric") ,
                 body = '
                   arma::mat v = Rcpp::as<arma::mat>(vs);
                   arma::mat op = v.i();
                   return Rcpp::List::create (Rcpp::Named ("outer") = op );',
                 plugin = "RcppArmadillo")

g2 <- cxxfunction(signature (vs = "numeric") ,
                  body = '
                   arma::mat v = Rcpp::as<arma::mat>(vs);
                   arma::mat op = v * v.t();
                   return Rcpp::List::create (Rcpp::Named ("outer") = op );',
                  plugin = "RcppArmadillo")

# start_time <- Sys.time()
# a<-g1(mat)
# end_time <- Sys.time()
# timeB1<-end_time-start_time
# 
# start_time <- Sys.time()
# b<-g2(mat)
# end_time <- Sys.time()
# timeB2<-end_time-start_time
# 
# c(timeB1,timeB2)

### Benchmark for inverse

benchmark("sparse" = {
  data <- rnorm(1e6)
  zero_index <- sample(1e6)[1:5e5]
  data[zero_index] <- 0
  mat <- matrix(data, ncol=1000)
  mat_sparse <- Matrix(mat, sparse=TRUE)
  res <- solve(mat_sparse)
},
"R base" = {
  data <- rnorm(1e6)
  zero_index <- sample(1e6)[1:5e5]
  data[zero_index] <- 0
  mat <- matrix(data, ncol=1000)
  res <- solve(mat)
},
"armadillo" = {
  data <- rnorm(1e6)
  zero_index <- sample(1e6)[1:5e5]
  data[zero_index] <- 0
  mat <- matrix(data, ncol=1000)
  res <- g1(mat)
},
replications = 100,
columns = c("test", "replications", "elapsed",
            "relative", "user.self", "sys.self"))


# test         replications elapsed relative user.self sys.self
# 3 armadillo          100  91.579    1.000    86.752    1.088
# 2    R base          100 124.615    1.361   113.763    1.251
# 1    sparse          100 242.292    2.646   202.983    9.732

### Benchmark for matrix multiplication

benchmark("sparse" = {
  data <- rnorm(1e6)
  zero_index <- sample(1e6)[1:5e5]
  data[zero_index] <- 0
  mat <- matrix(data, ncol=1000)
  mat_sparse <- Matrix(mat, sparse=TRUE)
  res <- mat_sparse %*% mat_sparse
},
"R base" = {
  data <- rnorm(1e6)
  zero_index <- sample(1e6)[1:5e5]
  data[zero_index] <- 0
  mat <- matrix(data, ncol=1000)
  res <- mat %*% mat
},
"armadillo" = {
  data <- rnorm(1e6)
  zero_index <- sample(1e6)[1:5e5]
  data[zero_index] <- 0
  mat <- matrix(data, ncol=1000)
  res <- g2(mat)
},
replications = 100,
columns = c("test", "replications", "elapsed",
            "relative", "user.self", "sys.self"))

# test         replications elapsed relative user.self sys.self
# 3 armadillo          100  40.447    1.000    37.899    1.021
# 2    R base          100  93.846    2.320    91.045    1.055
# 1    sparse          100  95.268    2.355    81.732    2.491