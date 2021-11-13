
######### individual-level randomization ########

#################### HTE ############################
L1_HTE_nc <- function(eff, alpha0, sigma2x=1, sigma2y=1, beta=0.2, alpha=0.05, m, ns, pi=0.5){
  
  lambda1 <- 1-alpha0
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1)/(sigma2w*sigma2x*ns*m)
  nc <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*sigma24/eff^2
  nc <- 2*(ceiling(nc/2))
  return (nc)
}


L1_HTE_power <- function(eff, alpha0, sigma2x=1, sigma2y=1, nc, alpha=0.05, m, ns, pi=0.5){
  
  lambda1 <- 1-alpha0
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1)/(sigma2w*sigma2x*ns*m)
  power <- pnorm( sqrt(nc*eff^2/sigma24)-qnorm(1-alpha/2) )
  return (power)
}
#######################################################




#################### ATE ############################

L1_ATE_nc <- function(eff, alpha0, sigma2y=1, beta=0.2, alpha=0.05, m, ns, pi=0.5){
  
  lambda1 <- 1-alpha0
  nc <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*sigma2y*lambda1/(pi*(1-pi)*ns*m*eff^2)
  nc <- 2*(ceiling(nc/2))
  return (nc)
}

L1_ATE_power <- function(eff, alpha0, sigma2y=1, nc, alpha=0.05, m, ns, pi=0.5){
  
  lambda1 <- 1-alpha0
  power <- pnorm( sqrt(nc*pi*(1-pi)*ns*m*eff^2/(sigma2y*lambda1))-qnorm(1-alpha/2) )
  return (power)
}
#####################################################




################# unadjusted ATE #####################

L1_ATE2_nc <- function(beta3, beta4, eff, alpha0, rho0, sigma2y=1, sigma2x=1, beta=0.2, alpha=0.05, m, ns, pi=0.5){
  
  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0
  
  lambda1 <- 1 - unadjusted_alpha0
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  nc <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*unadjusted_sigma2y*lambda1/(pi*(1-pi)*ns*m*eff^2)
  nc <- 2*(ceiling(nc/2))
  return (nc)
}

L1_ATE2_power <- function(beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y=1, sigma2x=1, nc, alpha=0.05, m, ns, pi=0.5){
  
  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0
  
  lambda1 <- 1 - unadjusted_alpha0
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  power <- pnorm( sqrt(nc*pi*(1-pi)*ns*m*eff^2/(unadjusted_sigma2y*lambda1))-qnorm(1-alpha/2) )
  return (power)
}

########################################################

