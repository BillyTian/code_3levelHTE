
######### subcluster-level randomization ########

#################### HTE ############################
L2_HTE_nc <- function(eff, rho0, alpha0, alpha1, sigma2x=1, sigma2y=1, beta=0.2, alpha=0.05, m, ns, pi=0.5){
  
  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(ns-1)*alpha1
  zeta1 <- 1-rho0
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1*lambda2*lambda3)/(sigma2w*sigma2x*(ns*m*lambda2*lambda3-ns*lambda3*(lambda2-lambda1)*(1+(m-1)*rho0)))
  nc <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*sigma24/eff^2
  nc <- 2*(ceiling(nc/2))
  return (nc)
}


L2_HTE_power <- function(eff, rho0, alpha0, alpha1, sigma2x=1, sigma2y=1, nc, alpha=0.05, m, ns, pi=0.5){

  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(ns-1)*alpha1
  zeta1 <- 1-rho0
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1*lambda2*lambda3)/(sigma2w*sigma2x*(ns*m*lambda2*lambda3-ns*lambda3*(lambda2-lambda1)*(1+(m-1)*rho0)))
  power <- pnorm( sqrt(nc*eff^2/sigma24)-qnorm(1-alpha/2) )
  return (power)
}
#######################################################




#################### ATE ############################
L2_ATE_nc <- function(eff, alpha0, alpha1, sigma2y=1, beta=0.2, alpha=0.05, m, ns, pi=0.5){
  
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  nc <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*sigma2y*lambda2/(pi*(1-pi)*ns*m*eff^2)
  nc <- 2*(ceiling(nc/2))
  return (nc)
}

L2_ATE_power <- function(eff, alpha0, alpha1, sigma2y=1, nc, alpha=0.05, m, ns, pi=0.5){

  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  power <- pnorm( sqrt(nc*pi*(1-pi)*ns*m*eff^2/(sigma2y*lambda2))-qnorm(1-alpha/2) )
  return (power)
}

#####################################################




################# unadjusted ATE #####################
L2_ATE2_nc <- function(beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y=1, sigma2x=1, beta=0.2, alpha=0.05, m, ns, pi=0.5){
  
  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0
  unadjusted_alpha1 <- w*alpha1 + (1-w)*rho1
  lambda2 <- 1 + (m-1)*unadjusted_alpha0 - m*unadjusted_alpha1
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  nc <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*unadjusted_sigma2y*lambda2/(pi*(1-pi)*ns*m*eff^2)
  nc <- 2*(ceiling(nc/2))
  return (nc)
}

L2_ATE2_power <- function(beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y=1, sigma2x=1, nc, alpha=0.05, m, ns, pi=0.5){
  
  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0
  unadjusted_alpha1 <- w*alpha1 + (1-w)*rho1
  lambda2 <- 1 + (m-1)*unadjusted_alpha0 - m*unadjusted_alpha1
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  power <- pnorm( sqrt(nc*pi*(1-pi)*ns*m*eff^2/(unadjusted_sigma2y*lambda2))-qnorm(1-alpha/2) )
  return (power)
}


########################################################