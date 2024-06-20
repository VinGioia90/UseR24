#################################################################
# Function for setting model formulas. Arguments:               #
# d: dimension of the outcome                                   #
# expl_mean - vector of strings for the mean vector predictors  #
# expl_Theta - vector of strings for the mean vector predictors #
#                                                               #
# Output: list containing the model formulas                    #
#################################################################

mformula <- function(d = 2, expl_mean, expl_Theta){ #, idxcov_int = rep(1, d * (d + 1)/2)){
  
  meanv_f <- paste0("y_", 1 : (d - 1), sep = "|", collapse = "") # First d - 1 outcome included in the model formula
  meanv_l <- paste0("y_", d, "~ s(", expl_mean[1], ") + s(", expl_mean[2],") + s(", expl_mean[3],")") # last part of the model formula (not beautiful:
  # only limited to this specific case)
  mean_foo <- formula(paste0(meanv_f, meanv_l)) # Complete mean model formula specification
  
  # Create the Theta labels (discriminating the intercepts from the covariate-dependent case)
  labelTh <- rep(0, d * (d + 1)/2)
  for(j in (d + 1) : (d + d * (d + 1)/2)) labelTh[j - d] <- SCM:::labTh(d, j)
  #labelTh_i <- labelTh[idxcov_int != 1]
  labelTh_x <- labelTh#[idxcov_int == 1]
  
  #lint <- sum(idxcov_int != 1)
  #lnoint <- sum(idxcov_int == 1)
  
  # f and l are meant for the first and lst part of the model formula, respectively
  # i and x are designed for capturing the possibility to consider intercepts and covariate-dependent formulas, respectively
  covi_foo_f <- covx_foo_f <- covx_foo_l <- covi_foo_l <- c()
  
  # Covariate-dependent case
  if( lnoint > 0){
    if(lnoint > 1) {
      covx_foo_f <- paste0(labelTh_x, sep = "|", collapse = "") #[1 : #(lnoint - 1)], sep = "|", collapse = "")
    }
    covx_foo_l <- paste0(labelTh_x[lnoint], "~ s(", expl_Theta[1], ") + s(", expl_Theta[2], ")") #not beautiful: only limited to this specific case
    covx_foo <- formula(paste0(covx_foo_f , covx_foo_l))
  }
  
  # Intercept case
  #if(lint > 0){
   # if(lint > 1){
  #    covi_foo_f <- paste0(labelTh_i[1 : (lint - 1)], sep = "|", collapse = "")
  #  }
  #  covi_foo_l <- paste0(labelTh_i[lint], "~ 1" )
  #  covi_foo <- formula(paste0(covi_foo_f , covi_foo_l))
  #}
  
  # Carrying out the final Theta model formula
  cov_foo <- unlist(covx_foo)
  #if(lnoint > 0 & lint == 0) cov_foo <- unlist(covx_foo)
  #if(lnoint == 0 & lint > 0) cov_foo <- unlist(covi_foo)
  #if(lnoint > 0 & lint > 0) cov_foo <- c(unlist(covx_foo), unlist(covi_foo))
  
  return(c(mean_foo, cov_foo))
}
