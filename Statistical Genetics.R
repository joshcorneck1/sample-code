################################################################
############ Code for a genetics simulation study ##############
################################################################
## Part (a)
# Generate the expression data
set.seed(01897838)
tic()
prop_diff <- c(0,0.1,0.3,0.6,0.9) # Vector of proportions
N <- 50 # Number of samples
m <- 20 # Number of genes in a pathway
K <- length(prop_diff) # Number of proportions
B <- 1000 # Number of permutations
L <- 30 # Number of runs
rej_f <- matrix(nrow=L,ncol=K) # Vector for Fisher rejections
rej_s <- matrix(nrow=L,ncol=K) # Vector for Stouffer rejections

for (l in 1:L){
  s <- runif(m,0.1,0.2) # Sample s
  mu <- rnorm(m,0,sqrt(s)) # Sample mu
  a <- rWishart(1,m,0.01*diag(1,m)) # Generate from Wishart to ensure positive semi-definite
  cov_mat <- sapply(1:m, function(i) a[i,,1]) # Convert into matrix
  gene_expression <- rmvnorm(N,mu,cov_mat) # Generate data

  ################################################################
  # Generate the response - one for each proportion, but each using
  # the same gene expression data
  Y <- list()
  for (k in 1:K){
    alpha <- 27 # Define alpha
    s <- sample(c(1:m),size = m*prop_diff[k]) # Sample associated genes
    beta <- rep(0,m) 
    beta[s] <- rnorm(length(s),0,sqrt(2)) # Sample beta for associated genes
    Y[[k]] <- alpha + gene_expression %*% beta + rnorm(N,1) # Generate response
  }
  ################################################################
  # Compute the test statistics
  for (k in 1:K){
  dat <- data.frame(y = Y[[k]],X = gene_expression) # Data frame
  # P-values for each gene
  p_vals <- sapply(1:m, function(i) summary(lm(y~dat[,i+1],data=dat))$coefficients[2,4]) 
  # Observed Fisher statistic
  stat_f_obs <- -2*sum(log(p_vals))
  # Observed Stouffer statistic
  stat_s_obs <- sum(qnorm(1 - p_vals))/sqrt(m)
  # Vectors to store permuted statistics
  stat_f_perm <- rep(0,B)
  stat_s_perm <- rep(0,B)
  for (b in 1:B){
    dat$y <- dat$y[sample(N,N)] # Permute responses
    # Compute new p-values
    p_vals <- sapply(1:m, function(i) summary(lm(y~dat[,i+1],data=dat))$coefficients[2,4])
    # Compute permuted test statistics
    stat_f_perm[b] <- -2*sum(log(p_vals))
    stat_s_perm[b] <- sum(qnorm(1 - p_vals))/sqrt(m)
  }
  # Store p-values for each test for given k
  p_f <- (sum(stat_f_perm > stat_f_obs) + 1)/(B+1)
  p_s <- (sum(stat_s_perm > stat_s_obs) + 1)/(B+1)
  
  # Store if H_0 is rejected for each test
  rej_f[l,k] <- ifelse(p_f<0.05/K,1,0)
  rej_s[l,k] <- ifelse(p_s<0.05/K,1,0)
  }
}
toc()

rej_f_store <- rej_f
rej_s_store <- rej_s

# Produce the plots
library(ggplot2)
dat <- data.frame(x=c(0,10,30,60,90),rej=c(colMeans(rej_f_store),colMeans(rej_s_store)),type=c(rep("Fisher",5),rep("Stouffer",5)))
gg <- ggplot(data=dat) + geom_point(aes(x=x, y=rej, color=type)) + geom_line(aes(x=x, y=rej, color=type))+ theme_bw() + theme(text = element_text(size = 20),
panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black")) +
xlab("Proportion of significant genes (%)") + scale_x_continuous(breaks=c(0,10,30,60,90)) +
  ylab("Power") + labs(color="Method")

################################################################
## Part (c)(i)
# Generate the expression data
set.seed(01897838)
tic()
prop_diff <- c(0.1,0.3,0.6,0.9) # Vector of proportions
N <- 50 # Number of samples
m_i <- c(10,30,40,50) # Number of genes in a pathway
K <- length(prop_diff) # Number of proportions
B <- 1000 # Number of permutations
L <- 20 # Number of runs
rej_f <- matrix(nrow=L,ncol=K) # Vector for Fisher rejections
rej_s <- matrix(nrow=L,ncol=K) # Vector for Stouffer rejections
rej_f_ii <- matrix(nrow=length(m_i),ncol=K)
rej_s_ii <- matrix(nrow=length(m_i),ncol=K)
r <- 0
for (m in m_i){
  r <- r + 1
  for (l in 1:L){
    s <- runif(m,0.1,0.2) # Sample s
    mu <- rnorm(m,0,sqrt(s)) # Sample mu
    a <- rWishart(1,m,0.01*diag(1,m)) # Generate from Wishart to ensure positive semi-definite
    cov_mat <- sapply(1:m, function(i) a[i,,1]) # Convert into matrix
    gene_expression <- rmvnorm(N,mu,cov_mat) # Generate data
    
    ################################################################
    # Generate the response - one for each proportion, but each using
    # the same gene expression data
    Y <- list()
    for (k in 1:K){
      alpha <- 27 # Sample alpha
      s <- sample(c(1:m),size = m*prop_diff[k]) # Sample associated genes
      beta <- rep(0,m) 
      beta[s] <- rnorm(length(s),0,sqrt(2)) # Sample beta for associated genes
      Y[[k]] <- alpha + gene_expression %*% beta + rnorm(N,1) # Generate response
    }
    ################################################################
    # Compute the test statistics
    for (k in 1:K){
      dat <- data.frame(y = Y[[k]],X = gene_expression) # Data frame
      # P-values for each gene
      p_vals <- sapply(1:m, function(i) summary(lm(y~dat[,i+1],data=dat))$coefficients[2,4]) 
      # Observed Fisher statistic
      stat_f_obs <- -2*sum(log(p_vals))
      # Observed Stouffer statistic
      stat_s_obs <- sum(qnorm(1 - p_vals))/sqrt(m)
      # Vectors to store permuted statistics
      stat_f_perm <- rep(0,B)
      stat_s_perm <- rep(0,B)
      for (b in 1:B){
        dat$y <- dat$y[sample(N,N)] # Permute responses
        # Compute new p-values
        p_vals <- sapply(1:m, function(i) summary(lm(y~dat[,i+1],data=dat))$coefficients[2,4])
        # Compute permuted test statistics
        stat_f_perm[b] <- -2*sum(log(p_vals))
        stat_s_perm[b] <- sum(qnorm(1 - p_vals))/sqrt(m)
      }
      # Store p-values for each test for given k
      p_f <- (sum(stat_f_perm > stat_f_obs) + 1)/(B+1)
      p_s <- (sum(stat_s_perm > stat_s_obs) + 1)/(B+1)
      
      # Store if H_0 is rejected for each test
      rej_f[l,k] <- ifelse(p_f<0.05/K,1,0)
      rej_s[l,k] <- ifelse(p_s<0.05/K,1,0)
    }
  }
  rej_f_ii[r,] <- colMeans(rej_f)
  rej_s_ii[r,] <- colMeans(rej_s)
}
toc()
 
power_N_f <- rbind(rej_f_i[1:3,],colMeans(rej_f_store),rej_f_i[4,])
power_N_s <- rbind(rej_s_i[1:3,],colMeans(rej_s_store),rej_s_i[4,])

power_m_f <-  rbind(rej_f_ii[1,],colMeans(rej_f_store)[-1],rej_f_ii[2:4,])
power_m_s <-  rbind(rej_s_ii[1,],colMeans(rej_s_store)[-1],rej_s_ii[2:4,])

# Produce plots
dat <- data.frame(x=c(10,20,30,40,50),rej=as.vector(power_m_s),size=c(rep("10",5),rep("30",5),rep("60",5),rep("90",5)))
gg <- ggplot(data=dat) + geom_point(aes(x=x, y=rej, color=size)) + geom_line(aes(x=x, y=rej, color=size))+ theme_bw() + theme(text = element_text(size = 20),
      panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black")) +
      xlab("Number of genes in pathway") + scale_x_continuous(breaks=c(10,20,30,40,50)) +
      ylab("Power") + labs(color="Proportion (%)")