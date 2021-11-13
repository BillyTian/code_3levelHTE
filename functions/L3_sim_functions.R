#Functions used for simulation corresponding to level-3 randomization

##########################################################################################
#
# Function to calculate sample size for testing HTE under cluster-level randomization
#
##########################################################################################
L3_HTE_n <- function(eff, rho0, rho1, alpha0, alpha1, sigma2x=1, sigma2y=1, beta=0.2, alpha=0.05, m, p, pi=0.5){
  # Argument:
  # eff: effect size
  # rho0: within-provider covariate ICC
  # rho1: between-provider covariate ICC
  # alpha0: within-provider outcome ICC
  # alpha1: between-provider outcome ICC
  # sigma2x: covariate variance
  # sigma2y: outcome variance
  # beta: type II error
  # alpha: type I error
  # m: number of patients under each provider
  # p: number of providers under each cluster
  # pi: proportion of treated
  
  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(p-1)*alpha1
  l1 <- 1-rho0
  l2 <- 1+(m-1)*rho0-m*rho1
  l3 <- 1+(m-1)*rho0+m*(p-1)*rho1
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1*lambda2*lambda3)/(sigma2w*sigma2x*(p*(m-1)*lambda2*lambda3*l1+(p-1)*lambda1*lambda3*l2+lambda1*lambda2*l3))
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*sigma24/eff^2
  n <- 2*(ceiling(n/2))
  return (n)
}

##########################################################################################
#
# Function to calculate power for testing HTE under cluster-level randomization
#
##########################################################################################
L3_HTE_power <- function(eff, rho0, rho1, alpha0, alpha1, sigma2x=1, sigma2y=1, n, alpha=0.05, m, p, pi=0.5){
  # Argument:
  # eff: effect size
  # rho0: within-provider covariate ICC
  # rho1: between-provider covariate ICC
  # alpha0: within-provider outcome ICC
  # alpha1: between-provider outcome ICC
  # sigma2x: covariate variance
  # sigma2y: outcome variance
  # n: number of clusters
  # alpha: type I error
  # m: number of patients under each provider
  # p: number of providers under each cluster
  # pi: proportion of treated
  
  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(p-1)*alpha1
  l1 <- 1-rho0
  l2 <- 1+(m-1)*rho0-m*rho1
  l3 <- 1+(m-1)*rho0+m*(p-1)*rho1
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1*lambda2*lambda3)/(sigma2w*sigma2x*(p*(m-1)*lambda2*lambda3*l1+(p-1)*lambda1*lambda3*l2+lambda1*lambda2*l3))
  power <- pnorm( sqrt(n*eff^2/sigma24)-qnorm(1-alpha/2) )
  return (power)
}

############################################################################################################################
#
# Function to calculate sample size and power for testing ATE (based on t-approximation) under cluster-level randomization
#
############################################################################################################################
L3_OTE_n_power <- function(eff, alpha0, alpha1, sigma2y=1, beta=0.2, alpha=0.05, m, p, pi=0.5){
  
  lambda3 <- 1+(m-1)*alpha0+m*(p-1)*alpha1
  nVar <- sigma2y*lambda3/(pi*(1-pi)*p*m)
  
  n <- 2
  power <- 0
  while (power < 1-beta){
    n <- n+2
    power <- pt(qt(1-alpha/2, n-2), n-2, ncp=eff/sqrt(nVar/n), lower.tail = F) + pt(qt(alpha/2, n-2), n-2, ncp=eff/sqrt(nVar/n))
  }
  return(c(n, power))
}

######################################################################################################################################
#
# Function to calculate sample size and power for testing unadjusted ATE (based on t-approximation) under cluster-level randomization
#
######################################################################################################################################
L3_OTE2_n_power <- function(beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y=1, sigma2x=1, beta=0.2, alpha=0.05, m, p, pi=0.5){
  
  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0
  unadjusted_alpha1 <- w*alpha1 + (1-w)*rho1
  lambda3 <- 1+(m-1)*unadjusted_alpha0+m*(p-1)*unadjusted_alpha1
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  
  nVar <- unadjusted_sigma2y*lambda3/(pi*(1-pi)*p*m)
  
  n <- 2
  power <- 0
  while (power < 1-beta){
    n <- n+2
    power <- pt(qt(1-alpha/2, n-2), n-2, ncp=eff/sqrt(nVar/n), lower.tail = F) + pt(qt(alpha/2, n-2), n-2, ncp=eff/sqrt(nVar/n))
  }
  return(c(n, power))
}

##########################################################################
#
# Function to generate correlated data under cluster-level randomization
#
##########################################################################
L3_datagen <- function(beta1, beta2, beta3, beta4, alpha0, alpha1, rho0, rho1, sigma2y=1, sigma2x=1, n, m, p){
  
  id <- seq(1, (n*p*m), 1)
  provider <-rep(1:(n*p), each=m)
  cluster <- rep(1:n, each=m*p)
  data <- data.frame(cbind(id, provider, cluster))
  
  #Randomize the cluster-level treatment W with equal allocation
  W.assign <- sample(1:n, n/2)
  data$W <- ifelse(data$cluster %in% W.assign, 1, 0)
  
  #Generate X (one continuous individual-level covariate)
  c <- rnorm(n*p*m, 0, sqrt(sigma2x*(1-rho0)))
  b <- rep(rnorm(n*p, 0, sqrt(sigma2x*(rho0-rho1))), each=m)
  a <- rep(rnorm(n, 0, sqrt(sigma2x*rho1)), each=m*p)
  data$X <- mu_X + a+b+c
  data$X_centered <- data$X - mean(data$X)
  
  #Generate Y under null and alternative
  epsilon <- rnorm(n*p*m, 0, sqrt(sigma2y*(1-alpha0)))
  u <- rep(rnorm(n*p, 0, sqrt(sigma2y*(alpha0-alpha1))), each=m)
  gamma <- rep(rnorm(n, 0, sqrt(sigma2y*alpha1)), each=m*p)
  
  data$Y <- beta1+beta2*data$W+beta3*data$X+beta4*data$W*data$X + gamma+u+epsilon
  
  return(data)
}