#Functions used for simulation corresponding to level-1 randomization


##########################################################################################
#
# Function to calculate sample size for testing HTE under individual-level randomization
#
##########################################################################################
L1_HTE_n <- function(eff, alpha0, sigma2x=1, sigma2y=1, beta=0.2, alpha=0.05, m, p, pi=0.5){
  # Argument:
  # eff: effect size
  # rho0: within-provider covariate ICC
  # alpha0: within-provider outcome ICC
  # alpha1: between-provider outcome ICC
  # sigma2x: covariate variance
  # sigma2y: covariate-adjusted outcome variance
  # beta: type II error
  # alpha: type I error
  # m: number of patients under each provider
  # p: number of providers under each cluster
  # pi: proportion of treated
  
  lambda1 <- 1-alpha0
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1)/(sigma2w*sigma2x*p*m)
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*sigma24/eff^2
  n <- 2*(ceiling(n/2))
  return (n)
}

######################################################################################
#
# Function to calculate power for testing HTE under individual-level randomization
#
######################################################################################
L1_HTE_power <- function(eff, alpha0, sigma2x=1, sigma2y=1, n, alpha=0.05, m, p, pi=0.5){
  # Argument:
  # eff: effect size
  # rho0: within-provider covariate ICC
  # alpha0: within-provider outcome ICC
  # alpha1: between-provider outcome ICC
  # sigma2x: covariate variance
  # sigma2y: covariate-adjusted outcome variance
  # n: number of clusters
  # alpha: type I error
  # m: number of patients under each provider
  # p: number of providers under each cluster
  # pi: proportion of treated
  
  lambda1 <- 1-alpha0
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1)/(sigma2w*sigma2x*p*m)
  power <- pnorm( sqrt(n*eff^2/sigma24)-qnorm(1-alpha/2) )
  return (power)
}

##########################################################################################
#
# Function to calculate sample size for testing ATE under individual-level randomization
#
##########################################################################################
L1_OTE_n <- function(eff, alpha0, sigma2y=1, beta=0.2, alpha=0.05, m, p, pi=0.5){
  # Argument:
  # eff: effect size (beta2)
  # alpha0: within-provider outcome ICC
  # alpha1: between-provider outcome ICC
  # sigma2y: outcome variance
  # beta: type II error
  # alpha: type I error
  # m: number of patients under each provider
  # p: number of providers under each cluster
  # pi: proportion of treated
  
  lambda1 <- 1-alpha0
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*sigma2y*lambda1/(pi*(1-pi)*p*m*eff^2)
  n <- 2*(ceiling(n/2))
  return (n)
}

##########################################################################################
#
# Function to calculate power for testing ATE under individual-level randomization
#
##########################################################################################
L1_OTE_power <- function(eff, alpha0, sigma2y=1, n, alpha=0.05, m, p, pi=0.5){
  # Argument:
  # eff: effect size
  # alpha0: within-provider outcome ICC
  # alpha1: between-provider outcome ICC
  # sigma2y: outcome variance
  # n: number of clusters
  # alpha: type I error
  # m: number of patients under each provider
  # p: number of providers under each cluster
  # pi: proportion of treated
  
  lambda1 <- 1-alpha0
  power <- pnorm( sqrt(n*pi*(1-pi)*p*m*eff^2/(sigma2y*lambda1))-qnorm(1-alpha/2) )
  return (power)
}

#####################################################################################################
#
# Function to calculate sample size for testing unadjusted ATE under individual-level randomization
#
#####################################################################################################
L1_OTE2_n <- function(beta3, beta4, eff, alpha0, rho0, sigma2y=1, sigma2x=1, beta=0.2, alpha=0.05, m, p, pi=0.5){
  # Argument:
  # beta3: true covariate effect
  # beta4: true effect of treatment effect heterogeneity
  # eff: effect size (beta2)
  # alpha0: within-provider outcome ICC
  # alpha1: between-provider outcome ICC
  # rho0: within-provider covariate ICC
  # rho1: between-provider covariate ICC
  # sigma2y: outcome variance
  # sigma2x: covariate variance
  # beta: type II error
  # alpha: type I error
  # m: number of patients under each provider
  # p: number of providers under each cluster
  # pi: proportion of treated
  
  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0
  
  lambda1 <- 1 - unadjusted_alpha0
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*unadjusted_sigma2y*lambda1/(pi*(1-pi)*p*m*eff^2)
  n <- 2*(ceiling(n/2))
  return (n)
}

#####################################################################################################
#
# Function to calculate power for testing unadjusted ATE under individual-level randomization
#
#####################################################################################################
L1_OTE2_power <- function(beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y=1, sigma2x=1, n, alpha=0.05, m, p, pi=0.5){
  # Argument:
  # beta3: true covariate effect
  # beta4: true effect of treatment effect heterogeneity
  # eff: effect size (beta2)
  # alpha0: within-provider outcome ICC
  # alpha1: between-provider outcome ICC
  # rho0: within-provider covariate ICC
  # rho1: between-provider covariate ICC
  # sigma2y: outcome variance
  # sigma2x: covariate variance
  # n: number of clusters
  # alpha: type I error
  # m: number of patients under each provider
  # p: number of providers under each cluster
  # pi: proportion of treated
  
  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0
  
  lambda1 <- 1 - unadjusted_alpha0
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  power <- pnorm( sqrt(n*pi*(1-pi)*p*m*eff^2/(unadjusted_sigma2y*lambda1))-qnorm(1-alpha/2) )
  return (power)
}

############################################################################
#
# Function to generate correlated data under individual-level randomization
#
############################################################################
L1_datagen <- function(beta1, beta2, beta3, beta4, rho0, rho1, alpha0, alpha1, sigma2x=1, sigma2y=1, n, m, p){
  
  id <- seq(1, (n*p*m), 1)
  provider <-rep(1:(n*p), each=m)
  cluster <- rep(1:n, each=m*p)
  data <- data.frame(cbind(id, provider, cluster))
  
  #Randomize the individual-level treatment W_ijk with equal allocation
  ind.trt <- NULL
  for (i in 1:(p*n)){
    #the provider sizes in simulation design are all even
    ind.trt_byprovider <- sample(c(rep(0, (m/2)), rep(1, (m/2))))
    ind.trt <- c(ind.trt, ind.trt_byprovider)
  }
  data$W <- ind.trt
  
  #Generate X (one continuous individual-level covariate)
  c <- rnorm(n*p*m, 0, sqrt(sigma2x*(1-rho0)))
  b <- rep(rnorm(n*p, 0, sqrt(sigma2x*(rho0-rho1))), each=m)
  a <- rep(rnorm(n, 0, sqrt(sigma2x*rho1)), each=m*p)
  data$X <- mu_X+a+b+c
  
  data$X_centered <- data$X-mean(data$X)
  
  #Generate Y under null and alternative
  epsilon <- rnorm(n*p*m, 0, sqrt(sigma2y*(1-alpha0)))
  u <- rep(rnorm(n*p, 0, sqrt(sigma2y*(alpha0-alpha1))), each=m)
  gamma <- rep(rnorm(n, 0, sqrt(sigma2y*alpha1)), each=m*p)
  
  data$Y <- beta1+beta2*data$W+beta3*data$X+beta4*data$W*data$X + gamma+u+epsilon
  
  return(data)
}