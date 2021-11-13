
######### Cluster-level randomization ########

#################### HTE ############################
L3_HTE_nc <- function(eff, rho0, rho1, alpha0, alpha1, sigma2x=1, sigma2y=1, beta=0.2, alpha=0.05, m, ns, pi=0.5){
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
  # ns: number of providers under each cluster
  # pi: proportion of treated
  
  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(ns-1)*alpha1
  l1 <- 1-rho0
  l2 <- 1+(m-1)*rho0-m*rho1
  l3 <- 1+(m-1)*rho0+m*(ns-1)*rho1
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1*lambda2*lambda3)/(sigma2w*sigma2x*(ns*(m-1)*lambda2*lambda3*l1+(ns-1)*lambda1*lambda3*l2+lambda1*lambda2*l3))
  nc <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*sigma24/eff^2
  nc <- 2*(ceiling(nc/2))
  return (nc)
}


L3_HTE_power <- function(eff, rho0, rho1, alpha0, alpha1, sigma2x=1, sigma2y=1, nc, alpha=0.05, m, ns, pi=0.5){
  
  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(ns-1)*alpha1
  l1 <- 1-rho0
  l2 <- 1+(m-1)*rho0-m*rho1
  l3 <- 1+(m-1)*rho0+m*(ns-1)*rho1
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1*lambda2*lambda3)/(sigma2w*sigma2x*(ns*(m-1)*lambda2*lambda3*l1+(ns-1)*lambda1*lambda3*l2+lambda1*lambda2*l3))
  power <- pnorm( sqrt(nc*eff^2/sigma24)-qnorm(1-alpha/2) )
  return (power)
}
#######################################################



#################### ATE (with t-approximation) ############################

L3_ATE_nc_power <- function(eff, alpha0, alpha1, sigma2y=1, beta=0.2, alpha=0.05, m, ns, pi=0.5){
  
  lambda3 <- 1+(m-1)*alpha0+m*(ns-1)*alpha1
  nVar <- sigma2y*lambda3/(pi*(1-pi)*ns*m)
  
  nc <- 2
  power <- 0
  while (power < 1-beta){
    nc <- nc+2
    power <- pt(qt(1-alpha/2, nc-2), nc-2, ncp=eff/sqrt(nVar/nc), lower.tail = F) + pt(qt(alpha/2, nc-2), nc-2, ncp=eff/sqrt(nVar/nc))
  }
  return(c(nc, power))
}
############################################################################



################# unadjusted ATE (with t-approximation) #####################

L3_ATE2_nc_power <- function(beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y=1, sigma2x=1, beta=0.2, alpha=0.05, m, ns, pi=0.5){
  
  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0
  unadjusted_alpha1 <- w*alpha1 + (1-w)*rho1
  lambda3 <- 1+(m-1)*unadjusted_alpha0+m*(ns-1)*unadjusted_alpha1
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  
  nVar <- unadjusted_sigma2y*lambda3/(pi*(1-pi)*ns*m)
  
  nc <- 2
  power <- 0
  while (power < 1-beta){
    nc <- nc+2
    power <- pt(qt(1-alpha/2, nc-2), nc-2, ncp=eff/sqrt(nVar/nc), lower.tail = F) + pt(qt(alpha/2, nc-2), nc-2, ncp=eff/sqrt(nVar/nc))
  }
  return(c(nc, power))
}
#############################################################################

