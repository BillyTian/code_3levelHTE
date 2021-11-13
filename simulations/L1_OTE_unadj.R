mainDir = '...' #change the directory
setwd(mainDir)

source("L1_sim_functions.R")

library(nlme)
delta_OTE <- 0.1 #beta2+mu_X*beta4
mu_X <- 1
beta3 <- 0.3
beta4 <- 0.05

m <- c(rep(20, 12), rep(50, 12))
p <- rep(c(rep(4, 6), rep(8, 6)), 2)
alpha0 <- rep(c(rep(0.015, 3), rep(0.1, 3)), 4)
alpha1 <- rep(c(rep(0.01, 3), rep(0.05, 3)), 4)
rho0 <- rep(c(0.15, 0.3, 0.5), 8)
rho1 <- rep(c(0.1, 0.15, 0.3), 8)
table <- cbind(m, p, alpha0, alpha1, rho0, rho1)

sigma2y <- 1
sigma2x <- 1
pi <- 0.5
Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
w <- sigma2y/(sigma2y+Q*sigma2x)
unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0
unadjusted_alpha1 <- w*alpha1 + (1-w)*rho1
table <- cbind(table, unadjusted_alpha0, unadjusted_alpha1)

n <- numeric(24)
pred.power <- numeric(24)
for (i in 1:nrow(table)){
  m.input <- as.numeric(table[i,][1])
  p.input <- as.numeric(table[i,][2])
  alpha0.input <- as.numeric(table[i,][3])
  #alpha1.input <- as.numeric(table[i,][4])
  rho0.input <- as.numeric(table[i,][5])
  #rho1.input <- as.numeric(table[i,][6])
  
  n[i] <- L1_OTE2_n(beta3=beta3, beta4=beta4, eff=delta_OTE, m=m.input, p=p.input, alpha0=alpha0.input, rho0=rho0.input)
  pred.power[i] <- L1_OTE2_power(beta3=beta3, beta4=beta4, n=n[i], eff=delta_OTE, m=m.input, p=p.input, alpha0=alpha0.input, rho0=rho0.input)
}

table <- cbind(table, n)


#function to compute empirical power or empirical type I error
empirical_OTE_unadjusted <- function(nullcase=F, parameter, nsims=5000){
  
  m <- as.numeric(parameter[1])
  p <- as.numeric(parameter[2])
  alpha0 <- as.numeric(parameter[3])
  alpha1 <- as.numeric(parameter[4])
  rho0 <- as.numeric(parameter[5])
  rho1 <- as.numeric(parameter[6])
  n <- as.numeric(parameter[9])
  
  if (nullcase==T){
    delta_OTE <- 0
  }
  
  beta1 <- 1
  beta2 <- delta_OTE-mu_X*beta4
  
  pvalue <- NULL
  count <- NULL
  est.beta <- NULL
  est.var <- NULL
  
  for (i in 1:nsims){
    set.seed(0807+2021*i)
    simdata <- L1_datagen(beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4, 
                          rho0=rho0, rho1=rho1, alpha0=alpha0, alpha1=alpha1, n=n, m=m, p=p)
      
    fit <- try(lme(Y ~ W, data=simdata, random= list(cluster = ~ 1, provider = ~ 1)), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1
    R <- c(0,1)
    beta <- fit$coef$fixed
    test.stat <- as.numeric((t(R)%*%beta)^2/(t(R)%*%vcov(fit)%*%R))
    pvalue[i] <- 1-pchisq(test.stat,1)
    est.beta[i] <- as.numeric(fit$coefficients$fixed[2])
    est.var[i] <- as.numeric(fit$varFix[2,2])
    
  }
  
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  mean.betahat <- mean(est.beta, na.rm=T)
  mean.varhat <- mean(est.var, na.rm=T)
  emp.var <- var(est.beta, na.rm=T)
  return(c(empirical, error.rate, mean.betahat, mean.varhat, emp.var))
}

#Compute empirical power (and error rate)
empirical.power <- NULL
for (i in 1:nrow(table)){
  empirical.power <- rbind(empirical.power, empirical_OTE_unadjusted(parameter=table[i,]))
}

#Compute empirical type I error (and error rate)
empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe <- rbind(empirical.tIe, empirical_OTE_unadjusted(nullcase=T, parameter=table[i,]))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))
names(result)[10:20] <- c("emp.tIe", "nconv.rate0", "mean.betahat0", "mean.varhat0", "emp.var0", "emp.power", "nconv.rate1", "mean.betahat1", "mean.varhat1", "emp.var1", "pred.power") 

write.table(result, paste0(mainDir, "/L1_OTE_unadj_eff0.1.csv"), sep =",", row.names=F)

