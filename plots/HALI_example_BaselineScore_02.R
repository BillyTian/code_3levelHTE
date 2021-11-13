###########################################
# Data example for 3-level HTE paper
###########################################
setwd("...") #change the directory
source("L3_plot_functions.R")
source("L2_plot_functions.R")

# (1) HALI: baseline spelling score
# (1.1) HTE
pi=0.5
m=25
ns=4
sigma2y=1
eff=0.12
sigma2x=1

alpha0=0.104 
alpha1=0.008

rho0=0.2 
rho1=0.1

L3_HTE_nc(eff=eff, rho0=rho0, rho1=rho1, alpha0=alpha0, alpha1=alpha1, 
          sigma2x=sigma2x, sigma2y=sigma2y, beta=0.2, alpha=0.05, m=m, ns=ns, pi=pi)
# nc=24

# what about power?
L3_HTE_power(eff=eff, rho0=rho0, rho1=rho1, alpha0=alpha0, alpha1=alpha1, 
             sigma2x=sigma2x, sigma2y=sigma2y, nc=24, alpha=0.05, m=m, ns=ns, pi=pi)
# 80.5% power

# (1.2) Adjusted ATE
L3_ATE_nc_power(eff=3*eff, alpha0=alpha0, alpha1=alpha1, 
                sigma2y=sigma2y, beta=0.2, alpha=0.05, m=m, ns=ns, pi=pi)

# (1.3) Unadjusted ATE
L3_ATE2_nc_power(beta3=eff*3, beta4=eff, eff=3*eff, alpha0=alpha0, alpha1=alpha1, rho0=rho0, rho1=rho1, 
                   sigma2y=sigma2y, sigma2x=sigma2x, beta=0.2, alpha=0.05, m=m, ns=ns, pi=pi)


pdf("HTE_data_1.pdf", width=12.5, height=12, paper="special")
par(mar=c(5.1, 6.1, 4.1, 2.1),mfrow=c(2,2))



### explore sensitivity with 24 clusters
### focus on HTE's and then comment on ATE's
x_seq<-seq(0,0.5,length=50)
y_seq<-seq(0,1,length=100)

mat_power<-outer(x_seq,y_seq,Vectorize(function(x,y)L3_HTE_power(eff=eff, rho0=x, rho1=(x*y), 
                                                                 alpha0=alpha0, alpha1=alpha1, 
                                                                 sigma2x=sigma2x, sigma2y=sigma2y, 
                                                                 nc=24, alpha=0.05, m=m, ns=ns, pi=pi)))
#levels=c(0.7,0.75,0.80,0.85,0.90,0.95),
contour(x_seq,y_seq,mat_power,labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(rho[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(rho[1],rho[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(a) HALI (Power with different covariate ICC)")),side=3,line=1,cex=1.5)



### explore sensitivity with 24 clusters
### focus on HTE's and then comment on ATE's
x_seq<-seq(0,0.2,length=50)
y_seq<-seq(0,1,length=100)

mat_power<-outer(x_seq,y_seq,Vectorize(function(x,y)L3_HTE_power(eff=eff, rho0=rho0, rho1=rho1, 
                                                                 alpha0=x, alpha1=(x*y), 
                                                                 sigma2x=sigma2x, sigma2y=sigma2y, 
                                                                 nc=24, alpha=0.05, m=m, ns=ns, pi=pi)))
#levels=c(0.7,0.75,0.80,0.85,0.90,0.95),
contour(x_seq,y_seq,mat_power,labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(alpha[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[1],alpha[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(b) HALI (Power with different outcome ICC)")),side=3,line=1,cex=1.5)


# (2) STRIDE: baseline spelling score
# (1.1) HTE
pi=0.5
m=63
ns=8
sigma2y=1
eff=0.2
sigma2x=0.16

alpha0=0.01
alpha1=0.005

rho0=0.1
# rho1=0.025

L2_HTE_nc(eff=eff, rho0=rho0, alpha0=alpha0, alpha1=alpha1, 
          sigma2x=sigma2x, sigma2y=sigma2y, beta=0.2, alpha=0.05, m=m, ns=ns, pi=pi)
# nc=10

# what about power?
L2_HTE_power(eff=eff, rho0=rho0, alpha0=alpha0, alpha1=alpha1, 
             sigma2x=sigma2x, sigma2y=sigma2y, nc=10, alpha=0.05, m=m, ns=ns, pi=pi)
# 80.4% power

# (1.2) Adjusted ATE
L2_ATE_nc(eff=(3/4)*eff, alpha0=alpha0, alpha1=alpha1, 
          sigma2y=sigma2y, beta=0.2, alpha=0.05, m=m, ns=ns, pi=pi)
L2_ATE_power(eff=(3/4)*eff, alpha0=alpha0, alpha1=alpha1, 
             sigma2y=sigma2y, alpha=0.05,nc=4,m=m, ns=ns, pi=pi)

# nc=4, 83.8%

### explore sensitivity with 80 clusters
### focus on HTE's and then comment on ATE's
x_seq<-seq(0,0.5,length=50)
y_seq<-seq(0,1,length=100)

mat_power<-outer(x_seq,y_seq,Vectorize(function(x,y)L2_HTE_power(eff=eff, rho0=x, 
                                                                 alpha0=alpha0, alpha1=alpha1, 
                                                                 sigma2x=sigma2x, sigma2y=sigma2y, 
                                                                 nc=10, alpha=0.05, m=m, ns=ns, pi=pi)))
#levels=c(0.7,0.75,0.80,0.85,0.90,0.95),
contour(x_seq,y_seq,mat_power,labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(rho[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(rho[1],rho[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(c) STRIDE (Power with different covariate ICC)")),side=3,line=1,cex=1.5)



### explore sensitivity with 80 clusters
### focus on HTE's and then comment on ATE's
x_seq<-seq(0,0.1,length=50)
y_seq<-seq(0,1,length=100)

mat_power<-outer(x_seq,y_seq,Vectorize(function(x,y)L2_HTE_power(eff=eff, rho0=rho0, 
                                                                 alpha0=x, alpha1=(x*y), 
                                                                 sigma2x=sigma2x, sigma2y=sigma2y, 
                                                                 nc=10, alpha=0.05, m=m, ns=ns, pi=pi)))
#levels=c(0.7,0.75,0.80,0.85,0.90,0.95),
contour(x_seq,y_seq,mat_power,labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(expression(alpha[0]),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[1],alpha[0])),side=2,line=3,cex=1.5,las=1)
mtext(expression(paste("(d) STRIDE (Power with different outcome ICC)")),side=3,line=1,cex=1.5)

dev.off()

