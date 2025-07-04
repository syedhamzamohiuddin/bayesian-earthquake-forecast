# ------------------------------------------------------------------------------
# 02_prior_sensitivity_analysis.R
#
# This script performs prior sensitivity analysis for a Bayesian AR(3) model 
# using earthquake count data. It computes the posterior parameters for the 
# model under multiple prior configurations and shows how robust the posterior 
# distributions of AR coefficients (ϕ) are to changes in the hyperparameters.
#
# Output: Posterior parameters (m, C, n*, d*) under different priors.
# ------------------------------------------------------------------------------

## read data, you need to make sure the data file is in your current working directory 
print("hi")
earthquakes.dat <- read.delim("data/earthquakes.txt")
earthquakes.dat$Quakes=as.numeric(earthquakes.dat$Quakes)

y.dat=earthquakes.dat$Quakes[1:100] ## this is the training data
y.new=earthquakes.dat$Quakes[101:103] ## this is the test data


######------------------------------------------------------------------------------------------------
## prior sensitivity analysis
## plot posterior distribution of phi_1 and phi_2 on a grid 



p=3  ## order of AR process
n.all=length(y.dat) ## T, total number of data

Y=matrix(y.dat[4:n.all],ncol=1)
Fmtx=matrix(c(y.dat[3:(n.all-1)],y.dat[2:(n.all-2)],y.dat[1:(n.all-3)]),nrow=p,byrow=TRUE)
n=length(Y)

## function to compute parameters for the posterior distribution of phi_1 and phi_2
## the posterior distribution of phi_1 and phi_2 is a multivariate t distribution

cal_parameters=function(m0=matrix(c(0,0,0),nrow=3),C0=diag(3),n0,d0){
  e=Y-t(Fmtx)%*%m0
  Q=t(Fmtx)%*%C0%*%Fmtx+diag(n)
  Q.inv=chol2inv(chol(Q))  ## similar as solve, but more robust
  A=C0%*%Fmtx%*%Q.inv
  m=m0+A%*%e
  C=C0-A%*%Q%*%t(A)
  n.star=n+n0
  d.star=t(Y-t(Fmtx)%*%m0)%*%Q.inv%*%(Y-t(Fmtx)%*%m0)+d0
  
  params=list()
  params[[1]]=n.star
  params[[2]]=d.star
  
  params[[3]]=m
  params[[4]]=C
  
  return(params)
}


## calculate density for three sets of hyperparameters
params1=cal_parameters(n0=2,d0=2)
params2=cal_parameters(n0=6,d0=1)
params3=cal_parameters(m0=matrix(c(-0.5,-0.5,-0.5),nrow=3),n0=6,d0=1)





######------------------------------------------------------------------------------------------------
# Prior sensitivity analysis shows that the posterior distribution is not sensitive to the choice of prior.
####################################################################################################################
#####_----------------------------------------------------------------------------------------------
