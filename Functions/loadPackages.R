# load packages 

remove.packages("mgcv")
if (!requireNamespace("mgcv", quietly = TRUE)) {
  print("Package not found, uninstallation was successful.")
  install.packages("/mgcv_developer/mgcv.tar.gz", repos = NULL, type="source")
  if(packageVersion("mgcv") != "9.0"){
    stop("Wrong version of mgcv!!")
  }
} else { 
  if(packageVersion("mgcv") != "9.0"){
    stop("Wrong version of mgcv!!")
  }
}

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


