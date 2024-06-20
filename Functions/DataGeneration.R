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
  
  
  ## Generating functions for the mean vectors and the Theta matrix
  f1m <- function(x1, m1) m1 * sin(pi * x1)
  f2m <- function(x2, m2) exp(m2 * x2)
  f3m <- function(x3, m3, m4, m5, m6) m3 * x3 ^ 11 * (m4 * (1 - x3)) ^ 6 + m5 * (m6 * x3) ^ 3 * (1 - x3) ^ 10
  fvc <- function(x1,x2,  vc0, vc1, vc2, vc3, vc4, vc5, vc6){
    vc0 + vc1 * sin(2 * pi * (x1 + vc2)) + vc3 * cos(2 * pi * (x1 + vc2)) + vc4 * sin(2 * pi * (x2 + vc5)) + vc6 * cos(2 * pi * (x2 + vc5))
  }
  
  # Setting the generating parameter values
  # For the mean vectors
  m1 <- runif(d, 1, 3)
  m2 <- runif(d, 1, 3)
  m3 <- runif(d, 0, 0.5)
  m4 <- runif(d, 9, 11)
  m5 <- runif(d, 9, 11)
  m6 <- runif(d, 9, 11)
  
  # For the entries of the (reparametrised) covariance matrix, that is for genereting the entries of Theta
  vc0 <- runif(d * (d + 1)/2, -0.25, 0.25)
  vc1 <- runif(d * (d + 1)/2, -0.5, 0.5)
  vc2 <- runif(d * (d + 1)/2, -0.5, 0.5)
  vc3 <- runif(d * (d + 1)/2, -0.25, 0.25)
  vc4 <- runif(d * (d + 1)/2, -1, 1)
  vc5 <- runif(d * (d + 1)/2, -1, 1)
  vc6 <- runif(d * (d + 1)/2, -0.5, 0.5)
  
  # Generation
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
    
    # Converting the generated Theta entries into the sensible lpi vector ordering
    lpi <- c(mu[[i]], diag(Theta[[i]]), Theta[[i]][upper.tri(Theta[[i]], diag = FALSE)])
    
    # MCD parametrisation: data generation
    LD <- SCM::internal()$mcd_LD(lpi, d)
    L <- LD
    diag(L) <- rep(1, d)
    D <- diag(LD)
    Sigma[[i]] <- crossprod(t(L) * D) # L%*%diag(D^2)%*%t(L)
    u <- mvnfast::rmvn(1, rep(0, d), diag(rep(1, d)))
    y[i,] <- t(mu[[i]] + L %*% t(u * D))
  }
  
  # Build the data data.frame
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


