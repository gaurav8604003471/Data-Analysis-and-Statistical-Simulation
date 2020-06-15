dmvnorm <- function (x, mean = rep(0, p), sigma = diag(p), log = FALSE) 
{
  if (is.vector(x)) 
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  if (!missing(mean)) {
    if (!is.null(dim(mean))) 
      dim(mean) <- NULL
    if (length(mean) != p) 
      stop("mean and sigma have non-conforming size")
  }
  if (!missing(sigma)) {
    if (p != ncol(sigma)) 
      stop("x and sigma have non-conforming size")
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                     check.attributes = FALSE)) 
      stop("sigma must be a symmetric matrix")
  }
  dec <- tryCatch(chol(sigma), error = function(e) e)
  if (inherits(dec, "error")) {
    x.is.mu <- colSums(t(x) != mean) == 0
    logretval <- rep.int(-Inf, nrow(x))
    logretval[x.is.mu] <- Inf
  }
  else {
    tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
    rss <- colSums(tmp^2)
    logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * 
                                                        pi) - 0.5 * rss
  }
  names(logretval) <- rownames(x)
  if (log) 
    logretval
  else exp(logretval)
}



#######################################################################################
set.seed(8)
X <- as.matrix(dat)
# Function calculates the EM estimates for any given C
GLMMforC <- function(X, C)
{
  n <- length(X[, 1])
  # tol <- 1e-10
  maxit <- 1e3
  corr <- 0
  ######## Starting values ###################
  ## pi are equally split over C
  
  pi.list <- rep(1/C, C)
  
  # The means for each C cannot be the same, since then the updates will be same. Hence adding random noise
  mu.mat <- matrix(c(0), nrow = 2, ncol = C)
  for( c in 1:C)
  { 
    mu.mat[1,c] <- mean(X[,1])  - ((C-1)*5) + (2*(c-1)*5)
    mu.mat[2,c] <- mean(X[,2]) - ((C-1)*1.25) + (2*(c-1)*1.25)
  }
  
  
  
  sigma.star <- matrix(c(70, 18, 18, 8), 2, 2)
  Sigma <- array(sigma.star, dim = c(2,2,1,C))
  
  iter <- 0
  diff <- 100
  theta <- c(pi.list, mu.mat, Sigma)
  
  current <- theta
  
  Ep <- matrix(0, nrow = n, ncol = C)  # 
  tol <- 1e-10
  
  while(diff > tol )
  {
    iter <- iter + 1
    
    for(c in 1:C)
    {
      Ep[ ,c] <- pi.list[c]*exp(log(dmvnorm(X, mu.mat[,c], Sigma[ , , , c])))
    }
    
    Ep <- Ep/rowSums(Ep)
    
    # M-step
    pi.list <- colMeans(Ep)
    
    for(c in 1:C)
    {
      mu.mat[ ,c] <- colSums(Ep[ ,c]%*%X) / sum(Ep[ ,c])
      
    }
    for(c in 1:C)
    {
      Sum = matrix(0, nrow = 2, ncol = 2)
      for( j in 1:n)
      {
        Sum = Sum + ((X[j,] - mu.mat[,c])%*%(t((X[j,] - mu.mat[,c])))*Ep[j,c])
      }
      Sigma[ , , , c] <-  Sum / sum(Ep[ ,c])   
      if(det(Sigma[ , , , c]) < 2.5e-8)
      {
        Sigma[ , , , c] <- Sigma[ , , , c] + diag(0.00001,2)
        corr <- corr +1
      }
    }
    
    theta <- c(pi.list, mu.mat, Sigma)
    diff <- max( abs(theta - current))
    
    current <- theta
  }
  
  return(list("mu.mat" = mu.mat, "Sigma" = Sigma, "pi.list" = pi.list, "Ep" = Ep))
}


###############################################
########### Question 1(a)
########### Estimate different parameter Values for C = 4

Estimate <- GLMMforC(X, C = 4)

Estimate$mu.mat 
Estimate$Sigma 
Estimate$pi.list 
labels <- apply(Estimate$Ep, 1,  which.max)
plot(X, col = labels)


##################################################
####### Question 1(b)
###### Using cross-validation with the negative log-likelihood as a loss function
###### choose the best model among the 5 models with C = (2, 3, 4, 5, 6). 
#### Functions needed for cross validation
# Finds the negative log likelihood


loglike <- function(X, C, mu.mat, Sigma, pi.list)
{
  n <-  length(X[ , 1])
  like <- 0
  for(i in 1:n)
  {  
    foo <- 0
    for(c in 1:C)
    {
      foo <- foo +  pi.list[c]*dmvnorm(t(X[i,]), mu.mat[,c], Sigma[ , , , c])
    }
    like <- like + log(foo)
  }    
  rtn <- -like
  return(rtn)
  
}


# My data set

n <- length(X[, 1])
potC <- 2:6

################################################
# 10-fold Cross-validation

permutation <- sample(1:n, replace = FALSE)
K <- 10

# Uneven folds, but that is ok
test.index <- split(permutation, rep(1:K, length = n, each = n/K))

potC <- 2:6
CV.errorLike <- numeric(length = length(potC))
CV.errorLike

for(c in 1:length(potC))
{
  foo3 <- 0
  foo5 <- 0
  for(k in 1:K)
  {
    X.train <- X[-test.index[[k]], ]
    
    X.test <- X[test.index[[k]], ]
    
    fitGMM <- GLMMforC(X = X.train, C = potC[c])
    
    foo3 <- foo3 + loglike(X = X.test, C = potC[c], fitGMM$mu.mat, fitGMM$Sigma, fitGMM$pi.list)
  }
  CV.errorLike[c] <- foo3/n
  
  
}

print(CV.errorLike) 


######### End of cross validations


fitC6 <- GLMMforC(X, C = 6)
fitC5 <- GLMMforC(X, C = 5)
fitC4 <- GLMMforC(X, C = 4)
fitC3 <- GLMMforC(X, C = 3)
fitC2 <- GLMMforC(X, C = 2)

fitC2$mu.mat 
fitC3$mu.mat
fitC4$mu.mat 
fitC5$mu.mat
fitC6$mu.mat

fitC2$Sigma 
fitC3$Sigma
fitC4$Sigma 
fitC5$Sigma
fitC6$Sigma


fitC2$pi.list 
fitC3$pi.list
fitC4$pi.list 
fitC5$pi.list
fitC6$pi.list

labels2 <- apply(fitC2$Ep, 1,  which.max)
labels3 <- apply(fitC3$Ep, 1,  which.max)
labels4 <- apply(fitC4$Ep, 1,  which.max)
labels5 <- apply(fitC5$Ep, 1,  which.max)
labels6 <- apply(fitC6$Ep, 1,  which.max)

plot(X, col = labels2)
plot(X, col = labels3)
plot(X, col = labels4)
plot(X, col = labels5)
plot(X, col = labels6)  




########################################################
######## Question 1(c)
######## Akaike Information Criterion (AIC)
######## choose the best model among the 5 models with C = (2, 3, 4, 5, 6) 
#### Functions needed for AIC
# AIC criteria


aic <- function(X, C, mu.mat, Sigma, pi.list)
{
  
  n <-  length(X[, 1])
  like <- 0
  foo <- 0
  for(i in 1:n)
  {  
    
    foo <- 0
    for(c in 1:C)
    {
      d <- length(X[ 1, ])
      Y <- matrix(0, nrow = 1, ncol = d)
      Y <- X[i,]
      foo <- foo +  pi.list[c]*dmvnorm(X[i,], mu.mat[,c], Sigma[ , , , c])
      
      # print( pi.list[c]*dmvnorm(X[i,], mu.mat[,c], Sigma[ , , , c]))
    }
    like <- like + log(foo)
  } 
  rtn <- - (2*like) + (2*((C - 1) + (3*C) + (2*C))) # No. of params  = 3*C -  1
  return(rtn)
  
}


foo4 <- 0
CV.errorAIC <- numeric(length = length(potC))
CV.errorAIC
for(c in 1:length(potC))
{
  fitGMM1 <- GLMMforC(X = X, C = potC[c])
  CV.errorAIC[c] <-  aic(X = X, C = potC[c], fitGMM1$mu.mat, fitGMM1$Sigma, fitGMM1$pi.list)
  
}
print(CV.errorAIC)



##########################################
######### Question 1(e)
######### Re-estimate the different parameter values for C = 5

Re_Estimate <- GLMMforC(X, C = 5)
Re_Estimate$mu.mat 
Re_Estimate$Sigma 
Re_Estimate$pi.list 


############################################
########## Question 1(f)
########## Scatterplot of data (X, Y ) for C = 5


labels1 <- apply(Re_Estimate$Ep, 1,  which.max)
plot(X, col = labels1)
