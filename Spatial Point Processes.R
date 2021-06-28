##############################################################
####### Code for MCMC on bivariate LGCP-Strauss model ########
##############################################################

# Load relevant packages
library(RandomFields)
library(spatstat)

# x and y coordinate vectors
x <- seq(0,1,length.out = 2^8)
y <- seq(0,1,length.out = 2^8)

# Define a function for computation of summary statistics
Tx <- function(X){
  X <- ppp(x=as.numeric(X[,1]),y=as.numeric(X[,2]),marks = as.factor(X[,3]))
  #1 
  n <- X$n
  L <- Lcross(X,correction = "iso")
  #2
  Lmax <- max(L$iso - L$r)
  #3
  Lmin <- min(L$iso - L$r)
  #4 
  Largmin <- L$r[which((L$iso - L$r) == max(L$iso - L$r))]
  #5
  L2 <- Lcross(X,correction = "iso",r=seq(0,0.2,length.out = m))
  Lspace <- L2$iso - L2$r
  #6
  q <- 5
  count <- quadratcount(X,nx=5,ny=5)
  Cmax <- max(count)/n
  #7
  Cmin <- min(count)/n
  #8
  countlist <- as.vector(count)/n
  Clogvar <- log(var(countlist))
  return(c(n,Lmax,Lmin,Largmin,Lspace,Cmax,Cmin,Clogvar))
}

# Number of pilot runs
kpilot <- 1000

# Matrix for parameters
par_samp <- matrix(nrow=kpilot,ncol=5)

# Matrix for summary statistics
Tsamp_plus <- matrix(nrow=kpilot,ncol=47)

# Maxmium number of iterations per run
maxit <- 20000

# Minimum number of observations to accept a run
m <- 40

# Do the pilot runs
for (k in 1:kpilot){
  # Counts for points of each mark
  N1 <- 0
  N2 <- 0
  # Run the algorithm (note requirement of N>= m)
  while(N1 < m || N2 < m){
    # Sample the parameters
    marks <- c("A","B")
    mu <- runif(1,3,6)
    sig2 <- runif(1,0.01,4)
    s <- runif(1,0.01,0.5)
    R <- runif(1,0,0.05)
    gam <- runif(1,0.01,1)
    
    # Save the parameters
    par_samp[k,] <- c(mu,sig2,s,R,gam)
    
    # Sample the field
    model <- RMexp(var=sig2,scale=s) + mu
    simu <- RFsimulate(model,x = seq(0,1,length.out = 2^8),y = seq(0,1,length.out = 2^8))
    z_vec <- c(exp(simu$variable1))
    z <- exp(RFspDataFrame2conventional(simu)$data)
    
    # Initialise with the point that maxmimises the field
    a <- which(z == max(z), arr.ind = TRUE)
    X <- matrix(c(x[a[2]],y[a[1]],sample(marks,1)),nrow=1,ncol=3)
    
    # Run the simulation for given values
    for (i in 1:maxit){
      mark <- sample(marks,1)
      p <- runif(1,0,1)
      if (p <= 0.5){
        # Birth proposal
        rx <- sample(c(1:length(x)),size=1)
        ry <- sample(c(1:length(y)),size=1)
        rm <- mark
        # Sampled point (and mark)
        u <- c(x[rx],y[ry],rm)
        # Points with a different mark
        diff_mark <- which(X[,3] != rm)
        # Compute number of near-points
        lR <- sum(unlist(sapply(1:length(diff_mark), function(j) dist(rbind(X[diff_mark[j],1:2],u[1:2])) < R)))
        if (is.na(lR)){
          lR <- 0
        }
        # Acceptance ratio
        A <- z[rx,ry]*gam^lR/(2*(nrow(X)+1))
        if (runif(1,0,1) < A){
          # Add the point
          X <- rbind(X,u)
        }
      }
      if (p > 0.5){
        # Death proposal
        r <- sample(1:nrow(X),1)
        u <- X[r,] 
        rx <- which(x==u[1])
        ry <- which(y==u[2])
        diff_mark <- which(X[,3] != u[3])
        lR <- sum(unlist(sapply(1:length(diff_mark), function(j) dist(rbind(X[diff_mark[j],1:2],u[1:2])) < R)))
        if (is.na(lR)){
          lR <- 0
        }
        # Acceptance ratio
        A <- 2*nrow(X)*(z[rx,ry]*gam^lR)^(-1)
        if (runif(1,0,1) < A && nrow(X) > 1){
          # Remove the point
          X <- matrix(X[-r,],ncol=3)
        }
      } 
    }
    N1 <- nrow(X[which(X[,3]=="A"),])
    N2 <- nrow(X[which(X[,3]=="B"),])
  }
  # Compute and save summary statistics
  Tsamp_plus[k,] <- Tx(X) 
  print(k)
}

# Save the data
save(Tsamp_plus,file= "summary_stats.RData")
save(par_samp,file = "params.RData")