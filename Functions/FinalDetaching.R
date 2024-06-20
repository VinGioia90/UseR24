remove.packages("mgcv")
if (!requireNamespace("mgcv", quietly = TRUE)) {
  print("Package not found, uninstallation was successful.")
  install.packages("mgcv", repos = NULL, type="source")
  if(packageVersion("mgcv") == "9.0"){
    stop("This is not the version on the CRAN!! Please proceed with manual disinstalling and installation")
  }
} else { 
  if(packageVersion("mgcv") == "9.0"){
    stop("This is not the version on the CRAN!! Please proceed with manual disinstalling and installation")
  }
}