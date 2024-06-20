# load packages 
library(mgcv)
library(SCM)
library(microbenchmark)

library(Rcpp)
sourceCpp(code = "
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
void over_writediag(arma::mat& E, const arma::vec& rho, const arma::vec& ind, const int k_sp) {
  int kk = k_sp - 1;
  int dd;
  for (int i = 0; i < ind.n_elem; ++i) {
    dd = ind(i) - 1;
    E(dd, dd) = exp(rho(kk) * 0.5);
  }
}
")

#################################################################
# Function for setting model formulas. Arguments:               #
# d: dimension of the outcome                                   #
# expl_mean - vector of strings for the mean vector predictors  #
# expl_Theta - vector of strings for the mean vector predictors #
#                                                               #
# Output: list containing the model formulas                    #
#################################################################

mformula <- function(d = 2, expl_mean, expl_Theta){ 
  
  meanv_f <- paste0("y_", 1 : (d - 1), sep = "|", collapse = "") 
  meanv_l <- paste0("y_", d, "~ s(", expl_mean[1], ") + s(", expl_mean[2],") + s(", expl_mean[3],")") 
  
  mean_foo <- formula(paste0(meanv_f, meanv_l)) 
  
  labelTh <- rep(0, d * (d + 1)/2)
  for(j in (d + 1) : (d + d * (d + 1)/2)) labelTh[j - d] <- SCM:::labTh(d, j)
  labelTh_x <- labelTh
  covi_foo_f <- covx_foo_f <- covx_foo_l <- covi_foo_l <- c()
  covx_foo_f <- paste0(labelTh_x[1:(length(labelTh_x) - 1) ], sep = "|", collapse = "") 
  covx_foo_l <- paste0(labelTh_x[length(labelTh_x) ], "~ s(", expl_Theta[1], ") + s(", expl_Theta[2], ")") 
  covx_foo <- formula(paste0(covx_foo_f , covx_foo_l))
  cov_foo <- unlist(covx_foo)
  
  return(c(mean_foo, cov_foo))
}

###################################################################
# Function for data generation. Arguments:                        #
#  d: dimension of the outcome                                    #
#  n: number of observations                                      #
#  seed: seed for generations                                     # 
#                                                                 #
# Output: list of                                                 #
#  data: (y,x) (matrix)                                           #
#  mu: simulated mu vectors (list of vectors)                     #
#  Sigma: simulated vcov (list of matrix)                         #
#  Theta: simulated Theta (list of matrix)                        #
###################################################################

datagen <- function(d = 2, n = 1000, seed = 13){
  
  set.seed(seed)
  y <- matrix(0, n, d)
  x <- list()
  x$x1 <- runif(n)
  x$x2 <- runif(n)
  x$x3 <- runif(n)
  
  f1m <- function(x1, m1) m1 * sin(pi * x1)
  f2m <- function(x2, m2) exp(m2 * x2)
  f3m <- function(x3, m3, m4, m5, m6) m3 * x3 ^ 11 * (m4 * (1 - x3)) ^ 6 + m5 * (m6 * x3) ^ 3 * (1 - x3) ^ 10
  fvc <- function(x1,x2,  vc0, vc1, vc2, vc3, vc4, vc5, vc6){
    vc0 + vc1 * sin(2 * pi * (x1 + vc2)) + vc3 * cos(2 * pi * (x1 + vc2)) + vc4 * sin(2 * pi * (x2 + vc5)) + vc6 * cos(2 * pi * (x2 + vc5))
  }
  
  m1 <- runif(d, 1, 3)
  m2 <- runif(d, 1, 3)
  m3 <- runif(d, 0, 0.5)
  m4 <- runif(d, 9, 11)
  m5 <- runif(d, 9, 11)
  m6 <- runif(d, 9, 11)
  
  vc0 <- runif(d * (d + 1)/2, -0.25, 0.25)
  vc1 <- runif(d * (d + 1)/2, -0.5, 0.5)
  vc2 <- runif(d * (d + 1)/2, -0.5, 0.5)
  vc3 <- runif(d * (d + 1)/2, -0.25, 0.25)
  vc4 <- runif(d * (d + 1)/2, -1, 1)
  vc5 <- runif(d * (d + 1)/2, -1, 1)
  vc6 <- runif(d * (d + 1)/2, -0.5, 0.5)
  
  mu <- list()
  Sigma <- list()
  Theta <- list()
  
  for(i in 1 : n){
    mu[[i]] <- rep(0, d)
    for(j in 1 : d) mu[[i]][j] <- f1m(x$x1[i], m1[j]) + f2m(x$x2[i], m2[j]) +  f3m(x$x3[i], m3[j], m4[j], m5[j], m6[j])
    
    Theta[[i]] <- matrix(0, d, d)
    count <- d + 1
    for(j in 1 : d){
      Theta[[i]][j, j] <- fvc(x$x1[i], x$x2[i], vc0[j], vc1[j], vc2[j], vc3[j], vc4[j], vc5[j], vc6[j])
      if(j > 1){
        for(k in 1 : (j - 1)){
          Theta[[i]][j, k] <- Theta[[i]][k, j] <- fvc(x$x1[i], x$x2[i], vc0[count], vc1[count], vc2[count], vc3[count], vc4[count], vc5[count], vc6[count])
          count <- count + 1
        }
      }
    }
    
    lpi <- c(mu[[i]], diag(Theta[[i]]), Theta[[i]][upper.tri(Theta[[i]], diag = FALSE)])
    
    LD <- SCM::internal()$mcd_LD(lpi, d)
    L <- LD
    diag(L) <- rep(1, d)
    D <- diag(LD)
    Sigma[[i]] <- crossprod(t(L) * D) 
    u <- mvnfast::rmvn(1, rep(0, d), diag(rep(1, d)))
    y[i,] <- t(mu[[i]] + L %*% t(u * D))
  }
  
  colnames(y) <-  paste0("y_", 1 : d)
  X <- matrix(0, n, length(x))
  for(j in 1 : length(x)) X[,j] <-  x[[j]]
  colnames(X) <- names(x)
  data <- cbind(y, X)
  
  return(list(data = data.frame(data),
              mu = mu,
              Theta = Theta,
              Sigma = Sigma))
}

#################################################################################
# Function for (simulating data and) fitting the models. Arguments:  a          #
# nobs: number of observations                                                  #
# dgrid: grid of orucome dimension values                                       #
# expl_mean: vector of strings for the mean vector predictors                   #
# expl_Theta: vector of strings for the mean vector predictors                  #
#                                                                               #
# Output: list containing                                                       #
# time: fitting times                                                           #
#################################################################################
sim_est_efs <- function(nobs, dgrid, expl_mean, expl_Theta){
  
  
  res_sim <- function(dgrid, nobs, expl_mean, expl_Theta){
    dss <- lapply(1 : length(dgrid), function(jj){
      out <- list()
      out$sim <- datagen(d = dgrid[jj], n = nobs, seed = 13 * dgrid[jj])
      out$foo <-  mformula(d = dgrid[jj], expl_mean = expl_mean, expl_Theta = expl_Theta)
      
      time <- microbenchmark(fit <- gam_scm(out$foo,
                                            family=mvn_scm(d = dgrid[jj]),
                                            optimizer= "efs",  data = as.data.frame(out$sim$data)), times = 1L)
      out$time <- time$time                 
      return(out)
    })
    return(dss)
  }  
  
  
  out_res_sim <- function(.x){
    out2 <- res_sim(dgrid = dgrid, nobs = nobs, expl_mean = expl_mean, expl_Theta = expl_Theta)
    return(list(gen = out2))
  }
  
  res <- list()
  nrun <- 1 
  res <- lapply(1 : nrun, out_res_sim)
  
  return(res)
}

################################################
# Extract the computational times (in minutes) #
################################################

fit_time <- function(obj, param = "mcd", dgrid = NULL, nrun = 1) {
  if(param == "mcd") param2 <- "MCD" 
  
  data_time <-  data.frame(unlist(lapply(1 : length(dgrid),
                                         function(z) unlist(lapply(1 : nrun, function(x) obj[[x]]$gen[[z]]$time / (1e9 * 60))))),
                           rep(dgrid, each = nrun))
  
  data_time <- cbind(data_time, rep(param2, each = length(dgrid)))
  colnames(data_time) <- c("time", "d", "Type")
  
  data <- aggregate(data_time$time, list(data_time$d), FUN = mean)
  data_min <- aggregate(data_time$time, list(data_time$d), FUN = min)
  data_max <- aggregate(data_time$time, list(data_time$d), FUN = max)
  
  out <- cbind(data, data_min[,2], data_max[,2])
  sum_res_time <- as.data.frame(out)
  sum_res_time <- cbind(sum_res_time, rep(param2, each = length(dgrid)))
  colnames(sum_res_time) <- c("d", "mean_time", "min_time", "max_time", "Type")
  
  return(list(sum_res = data.frame(sum_res_time),
              res = data.frame(data_time)))
}


