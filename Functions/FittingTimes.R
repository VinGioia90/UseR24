################################################
# Extract the computational times (in minutes) #
################################################

fit_time <- function(obj, param = "mcd", dgrid = NULL, nrun = 1) {
  if(param == "mcd") param2 <- "MCD" 
  
  data_time <-  data.frame(unlist(lapply(1 : length(dgrid),
                                         function(z) unlist(lapply(1 : nrun, function(x) obj[[x]]$gen[[z]]$time / (1e9 * 60))))),
                           rep(dgrid, each = nrun))
  
  # Padding with useful info
  data_time <- cbind(data_time, rep(param2, each = length(dgrid)))
  colnames(data_time) <- c("time", "d", "Type")
  
  # Summary
  data <- aggregate(data_time$time, list(data_time$d), FUN = mean)
  data_min <- aggregate(data_time$time, list(data_time$d), FUN = min)
  data_max <- aggregate(data_time$time, list(data_time$d), FUN = max)
  
  # Final objects (overall and summary)
  out <- cbind(data, data_min[,2], data_max[,2])
  sum_res_time <- as.data.frame(out)
  sum_res_time <- cbind(sum_res_time, rep(param2, each = length(dgrid)))
  colnames(sum_res_time) <- c("d", "mean_time", "min_time", "max_time", "Type")
  
  return(list(sum_res = data.frame(sum_res_time),
              res = data.frame(data_time),
              iter = data.frame(data_iter),
              sum_iter = data.frame(sum_res_iter)))
}