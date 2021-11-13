###############################################
# level 2 relationships
###############################################

setwd("...") #change the directory

library(ggplot2)
library(reshape2)
library(metR)


L2_HTE <- function(sigma2y=1, sigma2x=1, pi=0.5, m, alpha0, alpha1, rho0){
  
  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  
  sigma42_2 <- sigma2y/(pi*(1-pi)*sigma2x) * m/(m/lambda1 - (1+(m-1)*rho0)*(1/lambda1-1/lambda2))
  return(sigma42_2)
}

pdf("Remark1_VaryALPHA_lv2.pdf", width=12.5, height=12, paper="special")
par(mar=c(5.1, 6.1, 4.1, 2.1),mfrow=c(2,2))


### (2.1) fixing covariate-ICC, explore outcome-ICC
pi=0.5
m=25
ns=8
sigma2y=1
sigma2x=1

rho0=0.01

x_seq<-seq(0,0.5,length=100) # alpha0
y_seq<-seq(0,1,length=100) # alpha1/alpha0

mats<-outer(x_seq,y_seq,Vectorize(function(x,y)L2_HTE(rho0=rho0, alpha0=x, alpha1=(x*y),
                                                      sigma2x=sigma2x, sigma2y=sigma2y, m=m, pi=pi)))

contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(alpha[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[1],alpha[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(a) Variance (",rho[0]," = ", 0.01, ")")),side=3,line=1,cex=1.5)


### (2.2) fixing covariate-ICC, explore outcome-ICC
pi=0.5
m=25
ns=8
sigma2y=1
sigma2x=1

rho0=0.25

x_seq<-seq(0,0.5,length=100) # alpha0
y_seq<-seq(0,1,length=100) # alpha1/alpha0

mats<-outer(x_seq,y_seq,Vectorize(function(x,y)L2_HTE(rho0=rho0, alpha0=x, alpha1=(x*y),
                                                      sigma2x=sigma2x, sigma2y=sigma2y, m=m, pi=pi)))

contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(alpha[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[1],alpha[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(b) Variance (",rho[0]," = ", 0.25, ")")),side=3,line=1,cex=1.5)


### (2.3) fixing covariate-ICC, explore outcome-ICC
pi=0.5
m=25
ns=8
sigma2y=1
sigma2x=1

rho0=0.5

x_seq<-seq(0,0.5,length=100) # alpha0
y_seq<-seq(0,1,length=100) # alpha1/alpha0

mats<-outer(x_seq,y_seq,Vectorize(function(x,y)L2_HTE(rho0=rho0, alpha0=x, alpha1=(x*y),
                                                      sigma2x=sigma2x, sigma2y=sigma2y, m=m, pi=pi)))

contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(alpha[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[1],alpha[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(c) Variance (",rho[0]," = ", 0.5, ")")),side=3,line=1,cex=1.5)


### (2.4) fixing covariate-ICC, explore outcome-ICC
pi=0.5
m=25
ns=8
sigma2y=1
sigma2x=1

rho0=1

x_seq<-seq(0,0.5,length=100) # alpha0
y_seq<-seq(0,1,length=100) # alpha1/alpha0

mats<-outer(x_seq,y_seq,Vectorize(function(x,y)L2_HTE(rho0=rho0, alpha0=x, alpha1=(x*y),
                                                      sigma2x=sigma2x, sigma2y=sigma2y, m=m, pi=pi)))

contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(alpha[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[1],alpha[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(d) Variance (",rho[0]," = ", 1, ")")),side=3,line=1,cex=1.5)

dev.off()
