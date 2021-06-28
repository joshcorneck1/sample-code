###########################################################
######## Code for running FFBSm particle filter ###########
###########################################################

## Form of the model:

# X_n = \rho X_{n-1} + \sigma V_n 
# Y_n = \beta\exp(X_n/2)W_n

# with W_n,V_n \sim N(0,1) and X_0 \sim N(0, \sigma^2/(1-\rho^2))

# Define initial parameters for the model
tmax = 100; rho = 0.91; sigma = 1; beta = 0.5;

# Vectors of 0s for X and Y
X <- rep(0,tmax+1)
Y <- rep(0,tmax+1)

# Set seed and initialise
set.seed(01897838)
X[1] <- sigma/sqrt(1-rho^2)*rnorm(1,0,1)
Y[1] <- beta*exp(X[1]/2)*rnorm(1,0,1)

# Generate the "true" data
for (n in 2:(tmax + 1)){
  X[n] <- rho*X[n-1] + sigma*rnorm(1,0,1)
  Y[n] <- beta*exp(X[n]/2)*rnorm(1,0,1)
}

# Function for running FFBSm
FFBSm <- function(N,rho){
  set.seed(123)
  # Simulations for optimal proposal
  xu <- matrix(nrow = (tmax + 1),ncol = N)
  # Stored simulations and their weights for reversal
  xuFFB <- matrix(nrow = (tmax + 1),ncol = N)
  qFFB <- matrix(nrow = (tmax + 1),ncol = N)
  # Matrix for anscestors and ancestor weights
  xa <- matrix(nrow = (tmax + 1),ncol = N)
  Wa <- matrix(nrow = (tmax + 1),ncol = N)
  # Resampled simulations
  x <- matrix(nrow = (tmax + 1),ncol = N)
  # Unnormalised weights
  qq <- matrix(nrow = (tmax + 1),ncol = N)
  log_w <- matrix(nrow = (tmax + 1),ncol = N)
  # Normalised weights
  q <- matrix(nrow = (tmax + 1),ncol = N)
  # Anscestors index
  I <-  matrix(nrow = (tmax + 1),ncol = N)
  
  # Initialise (do the first step of algorithm)
  xu [1,] <- rnorm(N,0,1)
  qq[1,] <- exp(-.5*(beta*exp(xu[1,]/2))^(-2)*(Y[1])^2)/(sqrt(2*pi)*beta*exp(xu[1,]/2))
  # Normalise
  q[1,] <- qq[1,]/sum(qq[1,]) 
  # Store the approximations and weights for FFBSm
  # xuFFB[k,i] = X_{k-1}^i pre re-sampling
  xuFFB[1,] <- xu[1,]
  qFFB[1,] <- q[1,]
  # Do the resampling
  I[1,] <- sample(1:N,N,replace=TRUE,prob=q[1,])
  x[1,] <- xu[1,I[1,]]
  # Note xa[n,] are the ancestors of xuFFB[n+1,]. That is, we have
  # xa[k,i] = x_{k-1}^{a_k(i)}
  # Wa[k,i] is W_k^a_k(i), since we take W_0^a_0(i) = W_0^i
  xa[1,] <- x[1,]
  Wa[1,] <- qFFB[1,I[1,]]
  
  # Now do the remaining steps
  # Again, we use the log-sum-exp trick
  for (i in 2:(tmax + 1)){
    xu[i,] <- rho*x[i-1,] + sigma*rnorm(N,0,1)
    log_w[i,] <- -.5*(beta*exp(xu[i,]/2))^(-2)*(Y[i])^2 - log(sqrt(2*pi)*beta) - 0.5*xu[i,]
    
    # log exp sum trick
    offset <- max(log_w[i,])
    log_w[i,] <- log_w[i,] - offset
    
    qq[i,] <- exp(log_w[i,])
    q[i,] <- qq[i,]/sum(qq[i,])
    # Store the approximations and weights for FFBSm
    xuFFB[i,] <- xu[i,]
    qFFB[i,] <- q[i,]
    
    # Samples
    I[i,] <- sample(1:N,N,replace=TRUE,prob=q[i,])
    # Resampling
    x[i,] <- xu[i,I[i,]]
    xa[i,] <- x[i,]
    Wa[i,] <- qFFB[i,I[i,]]
    x[1:(i-1),] <- x[1:(i-1),I[i,]]
  }
  
  fx <- function(xa,xb){
    -0.5*(2*pi*sigma^2)^(-0.5)*exp(-0.5*(sigma)^(-2)*(xa - rho*xb)^2)
  }
  
  # Matrix for backwards weights. Note WnT[k,i] is W_{k-1|T}^i
  WnT <- matrix(nrow = (tmax + 1),ncol = N)
  # Matrix for adjusted weights. Note W_tilde[k,] is \tilde{W}_{k,k+1|T}
  W_tilde <- matrix(nrow = tmax,ncol = N)
  # Matrix for denominators. Note denom[n,i]=sum_{l=1}^N(W_0^lf_theta(X_{n+1}^i|X_n^l)).
  denom <-  matrix(nrow = tmax,ncol = N)
  
  # Initialise
  WnT[tmax + 1,] <- qFFB[tmax+1,]
  
  # Compute recursively backwards 
  for (k in tmax:1){
    denom[k,] <- sapply(1:N, function(j) sum(qFFB[k,]*fx(xuFFB[k+1,j],xuFFB[k,])))
    for (i in 1:N){
    num <-  WnT[k+1,]*fx(xuFFB[k+1,],xuFFB[k,i])
    total_sum <- sum(num/denom[k,])*qFFB[k,i]
    WnT[k,i] <- total_sum
    }
  }
  
  # Do the same for the tilde weights
  # The first entry is slightly different - I am taking the anscestor of X_0
  # to be X_0. W_tilde[tmax,i] is \tilde{W}_{T-1,T|T}^i
  # This should be \tilde{W}_{0,1|T}^i
  W_tilde[1,] <- qFFB[1,]*WnT[2,]*fx(xuFFB[2,],xa[1,])/denom[1,]
  
  for (n in 1:(tmax-1)){
      W_tilde[n+1,] <-  Wa[n+1,]*WnT[n+2,]*fx(xuFFB[n+2,],xa[n+1,])/denom[n+1,]
  }
  
  # Define function for s_k
  sk <- function(x1,x2){
    sigma^(-2)*x2*(x1 - rho*x2)
  }
  
  # Compute the approximation - first part is sum k=1 to k=T, second part is k=0 
  # Here I am just setting the ancestors of X_0^i = 0 and the weights to 1??
  hat_s <- sum(sapply(1:tmax, function(n) sum(W_tilde[n,]*sk(xuFFB[n+1,],xa[n,])))) + rho*sum(qFFB[1,]*(xuFFB[1,]^2  - 1/(1-rho^2)))
  
  expecF <- sapply(1:(tmax + 1), function(n) sum(WnT[n,]*xuFFB[n,]))
  expecSIR <- sapply(1:(tmax+1), function(i) mean(x[i,]))
  varF <- sapply(1:(tmax + 1), function(n) sum(WnT[n,]*xuFFB[n,]^2) -  sum(WnT[n,]*xuFFB[n,])^2)
  varSIR <- sapply(1:(tmax+1), function(i) mean(x[i,]^2) - (mean(x[i,]))^2)
  list_data <- list(expecF,expecSIR,varF,varSIR)
  return(hat_s)
}

## Plots of algorithm output

library(latex2exp)

# Expectation
plot(1:(tmax+1),X,type="l",ylab = TeX("E(X_{n}|Y_{0:T})"),xlab="n")
lines(1:(tmax+1),A[[1]],col="red")
lines(1:(tmax+1),A[[2]],col="blue")

# Variance
plot(1:(tmax+1),A[[4]],type="l",ylab = TeX("Var(X_{n}|Y_{0:T})"),xlab="n",col="blue")
lines(1:(tmax+1),A[[3]],col="red",type="l")

