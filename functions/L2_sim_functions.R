#Functions used for simulation corresponding to level-2 randomization

##########################################################################################
#
# Function to calculate sample size for testing HTE under subcluster-level randomization
#
##########################################################################################
L2_HTE_n <- function(eff, rho0, alpha0, alpha1, sigma2x=1, sigma2y=1, beta=0.2, alpha=0.05, m, p, pi=0.5){
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
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(p-1)*alpha1
  zeta1 <- 1-rho0
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1*lambda2*lambda3)/(sigma2w*sigma2x*(p*m*lambda2*lambda3-p*lambda3*(lambda2-lambda1)*(1+(m-1)*rho0)))
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*sigma24/eff^2
  n <- 2*(ceiling(n/2))
  return (n)
}

##########################################################################################
#
# Function to calculate power for testing HTE under subcluster-level randomization
#
##########################################################################################
L2_HTE_power <- function(eff, rho0, alpha0, alpha1, sigma2x=1, sigma2y=1, n, alpha=0.05, m, p, pi=0.5){
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
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(p-1)*alpha1
  zeta1 <- 1-rho0
  sigma2w <- pi*(1-pi)
  #sigma24 <- (sigma2y*lambda1*lambda2*lambda3)/(sigma2w*sigma2x*(p*(m-1)*lambda2*lambda3*zeta1+(p-1)*lambda1*lambda3*(1+(m-1)*rho0)+lambda1*lambda2*(1+(m-1)*rho0)))
  sigma24 <- (sigma2y*lambda1*lambda2*lambda3)/(sigma2w*sigma2x*(p*m*lambda2*lambda3-p*lambda3*(lambda2-lambda1)*(1+(m-1)*rho0)))
  power <- pnorm( sqrt(n*eff^2/sigma24)-qnorm(1-alpha/2) )
  return (power)
}

##########################################################################################
#
# Function to calculate sample size for testing ATE under subcluster-level randomization
#
##########################################################################################
L2_OTE_n <- function(eff, alpha0, alpha1, sigma2y=1, beta=0.2, alpha=0.05, m, p, pi=0.5){
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
  
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*sigma2y*lambda2/(pi*(1-pi)*p*m*eff^2)
  n <- 2*(ceiling(n/2))
  return (n)
}

##########################################################################################
#
# Function to calculate power for testing ATE under subcluster-level randomization
#
##########################################################################################
L2_OTE_power <- function(eff, alpha0, alpha1, sigma2y=1, n, alpha=0.05, m, p, pi=0.5){
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
  
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  power <- pnorm( sqrt(n*pi*(1-pi)*p*m*eff^2/(sigma2y*lambda2))-qnorm(1-alpha/2) )
  return (power)
}

#####################################################################################################
#
# Function to calculate sample size for testing unadjusted ATE under subcluster-level randomization
#
#####################################################################################################
L2_OTE2_n <- function(beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y=1, sigma2x=1, beta=0.2, alpha=0.05, m, p, pi=0.5){
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
  unadjusted_alpha1 <- w*alpha1 + (1-w)*rho1
  lambda2 <- 1 + (m-1)*unadjusted_alpha0 - m*unadjusted_alpha1
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*unadjusted_sigma2y*lambda2/(pi*(1-pi)*p*m*eff^2)
  n <- 2*(ceiling(n/2))
  return (n)
}

#####################################################################################################
#
# Function to calculate power for testing unadjusted ATE under subcluster-level randomization
#
#####################################################################################################
L2_OTE2_power <- function(beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y=1, sigma2x=1, n, alpha=0.05, m, p, pi=0.5){
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
  unadjusted_alpha1 <- w*alpha1 + (1-w)*rho1
  lambda2 <- 1 + (m-1)*unadjusted_alpha0 - m*unadjusted_alpha1
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  power <- pnorm( sqrt(n*pi*(1-pi)*p*m*eff^2/(unadjusted_sigma2y*lambda2))-qnorm(1-alpha/2) )
  return (power)
}

#############################################################################
#
# Function to generate correlated data under subcluster-level randomization
#
#############################################################################
L2_datagen <- function(beta1, beta2, beta3, beta4, alpha0, alpha1, rho0, rho1, sigma2y=1, sigma2x=1, n, m, p){
  
  id <- seq(1, (n*p*m), 1)
  provider <-rep(1:(n*p), each=m)
  cluster <- rep(1:n, each=m*p)
  data <- data.frame(cbind(id, provider, cluster))
  
  #Randomize the provider-level treatment W_ij with equal allocation
  provider.trt <- NULL
  for (i in 1:n){
    if (p%%2==0){ #if the number of providers is even
      provider.trt_assign <- sample(c(rep(0, (p/2)), rep(1, (p/2))))
    } else {
      provider.trt_assign <- sample(c(rep(0, floor(p/2)), rep(1, floor(p/2)), rbinom(1,1,0.5)))
    }
    provider.trt_bycluster <- rep(provider.trt_assign, each=m)
    provider.trt <- c(provider.trt, provider.trt_bycluster)
  }
  data$W <- provider.trt
  
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