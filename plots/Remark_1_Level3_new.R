###############################################
# level 3 relationships
###############################################

setwd("...") #change the directory

library(ggplot2)
library(reshape2)
library(metR)


L3_HTE <- function(sigma2y=1, sigma2x=1, pi=0.5, ns, m, alpha0, alpha1, rho0, rho1){
  
  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+(ns-1)*m*alpha1
  zeta1 <- 1-rho0
  zeta2 <- 1+(m-1)*rho0-m*rho1
  zeta3 <- 1+(m-1)*rho0+(ns-1)*m*rho1
  
  sigma42_3 <- sigma2y/(pi*(1-pi)*sigma2x) * (ns*m)/(zeta3/lambda3 + (ns-1)*zeta2/lambda2 + ns*(m-1)*zeta1/lambda1)
  return(sigma42_3)
}



pdf("Remark1_VaryRHO_lv3.pdf", width=12.5, height=12, paper="special")
par(mar=c(5.1, 6.1, 4.1, 2.1),mfrow=c(2,2))


### (1.1) fixing outcome-ICC, explore covariate-ICC
pi=0.5
m=25
ns=8
sigma2y=1
sigma2x=1

alpha0=0.05
alpha1=0.025

x_seq<-seq(0,1,length=100) # rho0
y_seq<-seq(0,1,length=100) # rho1/rho0

mats<-outer(x_seq,y_seq,Vectorize(function(x,y)L3_HTE(rho0=x, rho1=(x*y), alpha0=alpha0, alpha1=alpha1, 
                                                      sigma2x=sigma2x, sigma2y=sigma2y, m=m, ns=ns, pi=pi)))

contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(rho[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(rho[1],rho[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(a) Variance (",alpha[0]," = ", 0.05, ", ",alpha[1]," = ", 0.025, ")")),side=3,line=1,cex=1.5)



### (1.2) fixing outcome-ICC, explore covariate-ICC
pi=0.5
m=25
ns=8
sigma2y=1
sigma2x=1

alpha0=0.05
alpha1=0.05

x_seq<-seq(0,1,length=100) # rho0
y_seq<-seq(0,1,length=100) # rho1/rho0

mats<-outer(x_seq,y_seq,Vectorize(function(x,y)L3_HTE(rho0=x, rho1=(x*y), alpha0=alpha0, alpha1=alpha1, 
                                                      sigma2x=sigma2x, sigma2y=sigma2y, m=m, ns=ns, pi=pi)))

contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(rho[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(rho[1],rho[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(b) Variance (",alpha[0]," = ", 0.05, ", ",alpha[1]," = ", 0.05, ")")),side=3,line=1,cex=1.5)


### (1.3) fixing outcome-ICC, explore covariate-ICC
pi=0.5
m=25
ns=8
sigma2y=1
sigma2x=1

alpha0=0.1
alpha1=0.05

x_seq<-seq(0,1,length=100) # rho0
y_seq<-seq(0,1,length=100) # rho1/rho0

mats<-outer(x_seq,y_seq,Vectorize(function(x,y)L3_HTE(rho0=x, rho1=(x*y), alpha0=alpha0, alpha1=alpha1, 
                                                      sigma2x=sigma2x, sigma2y=sigma2y, m=m, ns=ns, pi=pi)))

contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(rho[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(rho[1],rho[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(c) Variance (",alpha[0]," = ", 0.1, ", ",alpha[1]," = ", 0.05, ")")),side=3,line=1,cex=1.5)

### (1.4) fixing outcome-ICC, explore covariate-ICC
pi=0.5
m=25
ns=8
sigma2y=1
sigma2x=1

alpha0=0.2
alpha1=0.05

x_seq<-seq(0,1,length=100) # rho0
y_seq<-seq(0,1,length=100) # rho1/rho0

mats<-outer(x_seq,y_seq,Vectorize(function(x,y)L3_HTE(rho0=x, rho1=(x*y), alpha0=alpha0, alpha1=alpha1, 
                                                      sigma2x=sigma2x, sigma2y=sigma2y, m=m, ns=ns, pi=pi)))

contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(rho[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(rho[1],rho[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(d) Variance (",alpha[0]," = ", 0.2, ", ",alpha[1]," = ", 0.05, ")")),side=3,line=1,cex=1.5)


dev.off()


pdf("Remark1_VaryALPHA_lv3.pdf", width=12.5, height=12, paper="special")
par(mar=c(5.1, 6.1, 4.1, 2.1),mfrow=c(2,2))


# ### (2.1) fixing covariate-ICC, explore outcome-ICC
# pi=0.5
# m=25
# ns=8
# sigma2y=1
# sigma2x=1
# 
# rho0=0.1
# rho1=0.1
# 
# x_seq<-seq(0,0.5,length=100) # alpha0
# y_seq<-seq(0,1,length=100) # alpha1/alpha0
# 
# mats<-ns*m*outer(x_seq,y_seq,Vectorize(function(x,y)L3_HTE(rho0=rho0, rho1=rho1, alpha0=x, alpha1=(x*y),
#                                                            sigma2x=sigma2x, sigma2y=sigma2y, m=m, ns=ns, pi=pi)))
# 
# contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
# mtext(expression(alpha[0]),side=1,line=3,cex=1.5)
# mtext(expression(frac(alpha[1],alpha[0])),side=2,line=3,cex=1.5,las=1)
# mtext(expression(paste("(a) Variance (",rho[0]," = ", 0.1, ", ",rho[1]," = ", 0.1, ")")),side=3,line=1,cex=1.5)



### (2.2) fixing covariate-ICC, explore outcome-ICC
pi=0.5
m=25
ns=8
sigma2y=1
sigma2x=1

rho0=0.25
rho1=0.1

x_seq<-seq(0,0.5,length=100) # alpha0
y_seq<-seq(0,1,length=100) # alpha1/alpha0

mats<-outer(x_seq,y_seq,Vectorize(function(x,y)L3_HTE(rho0=rho0, rho1=rho1, alpha0=x, alpha1=(x*y),
                                                      sigma2x=sigma2x, sigma2y=sigma2y, m=m, ns=ns, pi=pi)))

contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(alpha[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[1],alpha[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(a) Variance (",rho[0]," = ", 0.25, ", ",rho[1]," = ", 0.1, ")")),side=3,line=1,cex=1.5)


### (2.3) fixing covariate-ICC, explore outcome-ICC
pi=0.5
m=25
ns=8
sigma2y=1
sigma2x=1

rho0=0.5
rho1=0.5

x_seq<-seq(0,0.5,length=100) # alpha0
y_seq<-seq(0,1,length=100) # alpha1/alpha0

mats<-outer(x_seq,y_seq,Vectorize(function(x,y)L3_HTE(rho0=rho0, rho1=rho1, alpha0=x, alpha1=(x*y),
                                                      sigma2x=sigma2x, sigma2y=sigma2y, m=m, ns=ns, pi=pi)))

contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(alpha[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[1],alpha[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(b) Variance (",rho[0]," = ", 0.5, ", ",rho[1]," = ", 0.5, ")")),side=3,line=1,cex=1.5)


### (2.4) fixing covariate-ICC, explore outcome-ICC
pi=0.5
m=25
ns=8
sigma2y=1
sigma2x=1

rho0=0.75
rho1=0.5

x_seq<-seq(0,0.5,length=100) # alpha0
y_seq<-seq(0,1,length=100) # alpha1/alpha0

mats<-outer(x_seq,y_seq,Vectorize(function(x,y)L3_HTE(rho0=rho0, rho1=rho1, alpha0=x, alpha1=(x*y),
                                                      sigma2x=sigma2x, sigma2y=sigma2y, m=m, ns=ns, pi=pi)))

contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(alpha[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[1],alpha[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(c) Variance (",rho[0]," = ", 0.75, ", ",rho[1]," = ", 0.5, ")")),side=3,line=1,cex=1.5)


### (2.5) fixing covariate-ICC, explore outcome-ICC
pi=0.5
m=25
ns=8
sigma2y=1
sigma2x=1

rho0=1
rho1=0.5

x_seq<-seq(0,0.5,length=100) # alpha0
y_seq<-seq(0,1,length=100) # alpha1/alpha0

mats<-outer(x_seq,y_seq,Vectorize(function(x,y)L3_HTE(rho0=rho0, rho1=rho1, alpha0=x, alpha1=(x*y),
                                                      sigma2x=sigma2x, sigma2y=sigma2y, m=m, ns=ns, pi=pi)))

contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(alpha[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[1],alpha[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(d) Variance (",rho[0]," = ", 1, ", ",rho[1]," = ", 0.5, ")")),side=3,line=1,cex=1.5)



# ### (2.6) fixing covariate-ICC, explore outcome-ICC
# pi=0.5
# m=25
# ns=8
# sigma2y=1
# sigma2x=1
# 
# rho0=1
# rho1=1
# 
# x_seq<-seq(0,0.5,length=100) # alpha0
# y_seq<-seq(0,1,length=100) # alpha1/alpha0
# 
# mats<-ns*m*outer(x_seq,y_seq,Vectorize(function(x,y)L3_HTE(rho0=rho0, rho1=rho1, alpha0=x, alpha1=(x*y),
#                                                            sigma2x=sigma2x, sigma2y=sigma2y, m=m, ns=ns, pi=pi)))
# 
# contour(x_seq,y_seq,(mats),labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
# mtext(expression(alpha[0]),side=1,line=3,cex=1.5)
# mtext(expression(frac(alpha[1],alpha[0])),side=2,line=3,cex=1.5,las=1)
# mtext(expression(paste("(f) Variance (",rho[0]," = ", 1, ", ",rho[1]," = ", 1, ")")),side=3,line=1,cex=1.5)

dev.off()







