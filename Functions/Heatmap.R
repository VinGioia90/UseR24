heatmap_FitCov <- function(PredCov, d, range_var = NULL, range_corr = c(0, 1), label = "h", 
                           label_xaxis, label_yaxis, col_var, col_cor){
  
  diagDf <- data.frame(
    var1 = c(paste0(label, d : 1)),
    var2 = c(paste0(label, 1 : d)),
    Stdev =  sqrt(diag(PredCov)[d : 1])
  )
  
  sDf <- data.frame(X1 = rep(NA, d ^ 2), X2 = rep(NA, d ^ 2), Corr = rep(NA, d ^ 2)) 
  
  
  count <- 1
  for(i in 1 : d){
    for(j in 1 : d){
      if(i == j){
        sDf$X1[count] <- paste0(label, i)
        sDf$X2[count] <- paste0(label, (d - j) + 1)
        sDf$Corr[count] <- 0
      } 
      if(j < i){
        sDf$X1[count] <- paste0(label, i)
        sDf$X2[count] <- paste0(label, (d - j) + 1)
        sDf$Corr[count] <- NA
      }
      if(j > i){
        sDf$X1[count] <- paste0(label, i)
        sDf$X2[count] <- paste0(label, (d - j) + 1)
        sDf$Corr[count] <- PredCov[j, i]
      }
      count <- count + 1
    }
  }
  
  gg1 <- ggplot(sDf, aes(X1, X2)) + 
    geom_tile(aes(fill = Corr)) +
    scale_fill_gradientn(colors =  col_cor, limits = range_corr, na.value = "white") +
    new_scale_fill() +
    geom_tile(data = diagDf, aes(var1, var2, fill = Stdev)) +
    scale_fill_gradientn(colors = col_var) +
    theme(aspect.ratio = 1,
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          legend.key.width  = unit(1, "lines"),
          legend.key.height = unit(1, "lines"),
          legend.text=element_text(size=12),
          legend.title=element_text(size=15)) +
    scale_x_discrete(labels = label_xaxis)+
    scale_y_discrete(labels = label_yaxis)
  
  return(invisible(gg1))
  
}
