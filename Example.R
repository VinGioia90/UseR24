rm(list=ls())

library(rstudioapi)
root_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_dir)

######################################################################
# Code for reproducing the results provided at UseR!2024 Conference  #
######################################################################

######################################
#Hourly electricity demand modelling #
######################################

# load the data 
load("GEF14_d4.Rdata") 

# load the needed packages  
library(mgcViz)
library(mgcv)
# Uncomment the following line to install the SCM R package
# devtools::install_github("VinGioia90/SCM")
library(SCM)
library(ggplot2)
library(ggnewscale)
library(lubridate)

# Specify the dimension of the outcome 
d <- 4 

# Specify the mean model formula 
M_f <- list(l_h8 | l_h12 | l_h16 | l_h20  ~ dow + s(doy, bs = "cc"),
            l_h8 ~ l24_h8 + s(t95_h8), l_h12 ~ l24_h12 + s(t95_h12),
            l_h16 ~ l24_h16 + s(t95_h16), l_h20 ~ l24_h20 + s(t95_h20))

# Specify the mean + covariance matrix model formula 
MC_f <- c(M_f, list(
  Th_11 | Th_22 | Th_33 | Th_44 | Th_12 | Th_23 | Th_34  ~  dow  + s(doy, bs = "cc"),
  Th_11 ~ s(t95_h8), Th_22 ~ s(t95_h12), Th_33 ~ s(t95_h16), Th_44 ~ s(t95_h20)))

# Fit a static model (covariance matrix does not vary as function of the covariates)
fitS <- gam_scm(M_f, family = mvn_scm(d = 4), data = GEF14_d4)

# Fit a dynamic model (covariance matrix vary as function of the covariates)
fitD <- gam_scm(MC_f, family = mvn_scm(d = 4), data = GEF14_d4)

# Compare model performances 
AIC(fitS, fitD)

# residuals 
resS <- residuals(fitS, type = "response")
resD <- residuals(fitD, type = "response")

# Univariate residuals checks 
fitS <- getViz(fitS, nsim = 100, post = TRUE)
fitD <- getViz(fitD, nsim = 100, post = TRUE)

plS <- check1D(fitS, list("doy")) +
  l_gridCheck1D(gridFun = function(.x) {
    mcd(cov(.x))[1, 1]}, showReps = FALSE, n = 50)

plD <- check1D(fitD, list("doy")) +
  l_gridCheck1D(gridFun = function(.x) {
    mcd(cov(.x))[1, 1]}, showReps = FALSE, n = 50)

plS[[1]] 
plD[[1]] 

Sigma_pred_fitS <- predict(fitS, type = "response")
Sigma_predMat_fitS <-  Sigma_mat(Sigma_pred_fitS[,-c(1 : d)])[[1]]


# Analyse the ALEs of doy on the Var(y1) and the Corr(y1, y4)
plot(ALE(fitD, "doy", type = "response",
         oind = SCM:::sel_elem(4)(1, 1))) +  l_ciPoly() + l_fitLine() + 
  #annotate("text",  x = 230, y = 300, label = paste("Emp. Var. Residuals (h8) =", round(var(res[,1]),2)), vjust = 1, hjust = 1, cex = 5)
  annotate("text",  x = 230, y = 300, label = paste("Static Var. (h8) =", round(Sigma_predMat_fitS[1, 1],2)), vjust = 1, hjust = 1, cex = 5)
#mean(tapply(residuals(fitD, type="response")[,1], GEF14_d4$doy, var), na.rm=TRUE)

plot(ALE(fitD, "doy", type = "response",
         oind = SCM:::sel_elem(4)(4, 1))) +  l_ciPoly() + l_fitLine() + 
  #annotate("text",  x = 365/2, y = 0.25, label = "Cor(l_h8, l_h20)", vjust = 1, hjust = 1, cex = 5) + 
  annotate("text",  x = 225, y = 0.25, label = paste("Static Cor. (h8, h20) = ", round(Sigma_predMat_fitS[4, 1],2)), vjust = 1, hjust = 1, cex = 5) 

#for(i in 1:366){
# AAA[i] <- cor(res[GEF14_d4$doy==i,1],res[GEF14_d4$doy==i,4])
#}
#mean(AAA,na.rm =TRUE)

# Get the predicted values and extract the correlation matrix (padded with variances on the diagonal ) 

Sigma_pred_fitD <- predict(fitD, type = "response")
Sigma_predMat_fitD <- lapply(1 : nrow(GEF14_d4), function(x) Sigma_mat(Sigma_pred_fitD[,-c(1 : d)])[[x]])

# Select the days having minimum correlations between y1 and y4 and maximum variances of y1 
idx_min_CorrD <- which.min(unlist(lapply(1 : length(Sigma_predMat_fitD), function(x) unlist(Sigma_predMat_fitD[[x]][4, 1]))))
idx_max_VarD <- which.max(unlist(lapply(1 : length(Sigma_predMat_fitD), function(x) unlist(Sigma_predMat_fitD[[x]][1, 1]))))


# Heatmap
source("Functions/Heatmap.R")

# Settting graphical parameters
col_cor <- rev(colorspace::sequential_hcl(palette = "Blues 3", n = 100))
col_var <- rev(colorspace::sequential_hcl(palette = "Red", n = 10)[1 : 5])
label_xaxis <- c("h8","h12","h16","h20")
label_yaxis <- c("h20","h16","h12","h8")


days <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")


pl1 <- heatmap_FitCov(Sigma_predMat_fitD[[idx_min_CorrD]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h", 
                      label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)

date <- as.Date(paste(GEF14_d4[idx_min_CorrD, "year"], GEF14_d4[idx_min_CorrD, "doy"]), format = "%Y %j")
dow <- days[wday(date)]
pl1 +  annotate("text",  x=4, y = 4.25, label = paste(dow, date), vjust=1, hjust=1, cex = 7)   


pl2 <- heatmap_FitCov(Sigma_predMat_fitD[[idx_max_VarD]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h", 
                      label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
date <- as.Date(paste(GEF14_d4[idx_max_VarD,"year"], GEF14_d4[idx_max_VarD,"doy"]), format = "%Y %j")
dow <- days[wday(date)]
pl2 +  annotate("text",  x = 4, y = 4.25, label = paste(dow, date), vjust = 1, hjust = 1, cex = 7)   

##################################################
# Computational times under numerical experiment #
##################################################
rm(list=ls())

.rs.restartR()

library(rstudioapi)
root_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_dir)


source("Functions/loadPackages.R")
source("Functions/ModelFormula.R")
source("Functions/DataGeneration.R")
source("Functions/ModelFitting.R")
source("Functions/FittingTimes")

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

source("Functions/FinalDetaching.R")