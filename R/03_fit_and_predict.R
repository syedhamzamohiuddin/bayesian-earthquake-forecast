# ------------------------------------------------------------------------------
# 03_posterior_sampling_and_forecasting.R
#
# This script fits a Bayesian AR(3) model to earthquake count data using a 
# Normal-Inverse Gamma prior. It performs posterior sampling (ϕ and ν) and 
# generates 4-step-ahead forecasts. Outputs include posterior means and 
# predicted earthquake counts.
# ------------------------------------------------------------------------------

library(mvtnorm)

# 5000 Posterior samples of model parameters ϕ and ν

# Load data
earthquakes.dat <- read.delim("data/earthquakes.txt")
earthquakes.dat$Quakes <- as.numeric(earthquakes.dat$Quakes)
y.dat <- earthquakes.dat$Quakes[1:100]


p=3 ## order of AR process
n.all=length(y.dat) ## T, total number of data


Y=matrix(y.dat[4:n.all],ncol=1)
Fmtx=matrix(c(y.dat[3:(n.all-1)],y.dat[2:(n.all-2)],y.dat[1:(n.all-3)]),nrow=p,byrow=TRUE)
n=length(Y)


## posterior inference


## set the prior
m0=matrix(rep(0,p),ncol=1)
C0=10*diag(p)
n0=0.02
d0=0.02


## calculate parameters that will be reused in the loop
e=Y-t(Fmtx)%*%m0
Q=t(Fmtx)%*%C0%*%Fmtx+diag(n)
Q.inv=chol2inv(chol(Q))
A=C0%*%Fmtx%*%Q.inv
m=m0+A%*%e
C=C0-A%*%Q%*%t(A)
n.star=n+n0
d.star=t(Y-t(Fmtx)%*%m0)%*%Q.inv%*%(Y-t(Fmtx)%*%m0)+d0




n.sample=5000


## store posterior samples
nu.sample=rep(0,n.sample)
phi.sample=matrix(0,nrow=n.sample,ncol=p)




for (i in 1:n.sample) {
  set.seed(i)
  nu.new=1/rgamma(1,shape=n.star/2,rate=d.star/2)
  nu.sample[i]=nu.new
  phi.new=rmvnorm(1,mean=m,sigma=nu.new*C)
  phi.sample[i,]=phi.new
}
######------------------------------------------------------------------------------------------------
## After running the code, 5000 posterior samples of ϕ and ν are stored in phi.sample and 
## nu.sample.The h-step ahead prediction can be obtained by the following code:

## the prediction function


y_pred_h_step=function(s){
  h.step=4
  phi.cur=matrix(phi.sample[s,],ncol=1)
  nu.cur=nu.sample[s]
  y.cur=c(y.dat[n.all],y.dat[(n.all-1)],y.dat[(n.all-2)])
  y.pred=rep(0,h.step)
  for (i in 1:h.step) {
    mu.y=sum(y.cur*phi.cur)
    y.new=rnorm(1,mu.y,sqrt(nu.cur))
    y.pred[i]=y.new
    y.cur=c(y.new,y.cur)
    y.cur=y.cur[-length(y.cur)]
  }
  return(y.pred)
}
print(mean(phi.sample[,1]))
print(mean(phi.sample[,2]))

print(mean(nu.sample))

str(phi.sample)

o = sapply(1:5000, y_pred_h_step)  
print(str(o))


print(mean(o[1,]))
print(mean(o[2,]))
print(mean(o[3,]))
print(mean(o[4,]))
print(length(y.dat))
mean(o,)
