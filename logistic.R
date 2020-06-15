######################################
##### Bayesian logistic regression with MCMC
####### Question 2(c)
##### MCMC sampling with help of MetropolisHastings algorithm


set.seed(5)
Z = as.data.frame(dat)
D = as.matrix(dat)
y <- D[,1]
X <- D[, 2:6]

log.post <- function(y, X, beta, s2 = 100, mu = 0)
{
  n <- length(y)
  p <- dim(X)[2]
  Sum_term <- 0
  for( i in 1:n)
  {
    Sum_term <- Sum_term + log(1+exp(X[i,]%*%beta))
  }
  foo <- (-p/2) * log(s2) - ((t(beta-mu)%*%(beta-mu)/2)/s2) + sum(y*X%*%beta) - Sum_term 
  # print(foo)
  return(foo)
}

BLRmcmc <- function(y, X, T , mu = 0, start, s2 = 100, h1 , h2 , h3 , h4  , h5 )
{
  p <- dim(X)[2]
  
  beta <- matrix(0, nrow = T, ncol = p)
  
  
  beta[1,] <- start[1:p]   # Set starting values
  
  
  acc <- 0
  for(t in 2:T)
  {
    prop <- c(beta[t-1,]) + rnorm(p, sd = sqrt(c(h1,h2,h3,h4,h5)))   # Proposal
    
    
    # difference of log posteriors
    ratio <- log.post(y = y, X = X, beta = prop[1:p], s2 = s2, mu = mu) - log.post(y = y, X = X, beta = beta[t-1,], s2 = s2, mu = mu)
    U <- runif(1)
    
    if(U < exp(ratio))
    {
      acc <- acc + 1
      beta[t,] <- prop[1:p]
      
    } else{
      beta[t, ] <- beta[t-1, ]
      
    }
    
  }
  print(paste("Acceptance prob = ", acc/T))
  return(list("beta.est" = beta, acc.prop = acc/T))
}

lm.obj <- glm(formula = y ~ V3+V4+V5+V6 ,  Z, family = "binomial")

chain <- BLRmcmc(T = 1e5, y = y, X = X,  start = c(0, 0, 0, 0, 0),  h1 = 0.3, h2 = 0.1, h3 = 0.1, h4 = 0.1, h5 = 0.1)
post.mean <- c(colMeans(chain$beta.est))
post.quant1 <- c(quantile(chain$beta.est[,1], .05), quantile(chain$beta.est[,2], .05),quantile(chain$beta.est[,3], .05), quantile(chain$beta.est[,4], .05),quantile(chain$beta.est[,5], .05))
post.quant2 <- c(quantile(chain$beta.est[,1], .95), quantile(chain$beta.est[,2], .95),quantile(chain$beta.est[,3], .95), quantile(chain$beta.est[,4], .95),quantile(chain$beta.est[,5], .95))


# Estimates of intercept
c(post.quant1[1],post.mean[1], post.quant2[1])
lm.obj$coefficients[1]  # Compare with MLE

# Estimates of coefficient of X1
c(post.quant1[2], post.mean[2], post.quant2[2])
lm.obj$coefficients[2]  # Compare with MLE

# Estimates of coefficient of X2
c(post.quant1[3], post.mean[3], post.quant2[3])
lm.obj$coefficients[3]  # Compare with MLE

# Estimates of coefficient of X3
c(post.quant1[4], post.mean[4], post.quant2[4])
lm.obj$coefficients[4]  # Compare with MLE

# Estimates of coefficient of X4
c(post.quant1[5], post.mean[5], post.quant2[5])
lm.obj$coefficients[5]  # Compare with MLE




############################################
########## Question 2(d)

#Visualizing plots

plot.ts(chain$beta.est[,1])
plot.ts(chain$beta.est[,2])
plot.ts(chain$beta.est[,3])
plot.ts(chain$beta.est[,4])
plot.ts(chain$beta.est[,5])


acf(chain$beta.est[ ,1])
acf(chain$beta.est[ ,2])
acf(chain$beta.est[ ,3])
acf(chain$beta.est[ ,4])
acf(chain$beta.est[ ,5])


plot(density(chain$beta.est[,1]))
abline(v = c(post.mean[1], post.quant1[1], post.quant2[1]), col = c("red", "blue", "blue"))

plot(density(chain$beta.est[,2]))
abline(v = c(post.mean[2], post.quant1[2], post.quant2[2]), col = c("red", "blue", "blue"))

plot(density(chain$beta.est[,3]))
abline(v = c(post.mean[3], post.quant1[3], post.quant2[3]), col = c("red", "blue", "blue"))

plot(density(chain$beta.est[,4]))
abline(v = c(post.mean[4], post.quant1[4], post.quant2[4]), col = c("red", "blue", "blue"))

plot(density(chain$beta.est[,5]))
abline(v = c(post.mean[5], post.quant1[5], post.quant2[5]), col = c("red", "blue", "blue"))


######################################################
########### Question 2(e)
print(paste("Acceptance prob = ", chain$acc.prop))
print(c(colMeans(chain$beta.est)))
