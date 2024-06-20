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
sim_est_efs <- function(nobs, dgrid, expl_mean, expl_Theta)){

  
  res_sim <- function(dgrid, nobs, expl_mean, expl_Theta){
    dss <- lapply(1 : length(dgrid), function(jj){
      out <- list()
      # Generating data
      out$sim <- datagen(d = dgrid[jj], n = nobs, seed = 13 * dgrid[jj])
      # Setting the model formula
      out$foo <-  mformula(d = dgrid[jj], expl_mean = expl_mean, expl_Theta = expl_Theta)
      
      # Computational times 
      time <- microbenchmark(fit <- gam_scm(out$foo,
                                              family=mvn_scm(d = dgrid[jj], param = param2),
                                              optimizer= "efs",  data = as.data.frame(out$sim$data)), times = 1L)
      out$time <- time$time                 
      }
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
