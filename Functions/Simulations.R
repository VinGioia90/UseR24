##################################################
# Computational times under numerical experiment #
##################################################

#################################################################
# !!!!NOTE This step is quite fragile:                          #
# -) we are removing the CRAN version of mgcv                   #
# -) we install the developer version                           #
# -) we restore the CRAN version of mgcv                        #
#                                                               # 
# If you are not sure to control these steps, please trust that #
# the developer version is able to speed up the model fitting   #
# and don't run the following lines of code (you can see the    #
# results at line 51)                                           #
#################################################################
rm(list=ls())

# Removing mgcv R package
remove.packages("mgcv")

# Restart R session
.rs.restartR()

# Get wd
library(rstudioapi)
root_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_dir)

# Install developer version
devtools::install(paste0(root_dir, "/mgcv_developer/mgcv"))

# Check
if(packageVersion("mgcv") == 9.0){
  print("!!! OK !!! you can go to the next step")
} else {
  print("!!! NO !!! If you run the following code your computational times are grater than those reported")
}

# Load all the needed functions and packages
source("Utilities.R")

# Running the experiment
nrun <- 1
dgrid <- c(2,5,10,15,20)
nobs <- 10000

sim_mcdG_mcdF <- sim_est_efs(nobs, dgrid, expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"))

TIMES <- fit_time(sim_mcdG_mcdF, param = "mcd", dgrid = dgrid, nrun=1)

TIMES$sum_res  
#   d   mean_time    min_time    max_time Type
#1  2  0.09269453  0.09269453  0.09269453  MCD
#2  5  0.56419284  0.56419284  0.56419284  MCD
#3 10  3.53005337  3.53005337  3.53005337  MCD
#4 15 13.42818425 13.42818425 13.42818425  MCD
#5 20 46.18585052 46.18585052 46.18585052  MCD


remove.packages("mgcv")
install.packages("mgcv") #!!! Digit yes !!!

if(packageVersion("mgcv") == 9.0){
  print("!!! NO !!! Your version is the developer version. Find a way to get the CRAN version")
} else {
  print("OK. The process ends succesfully")
}



