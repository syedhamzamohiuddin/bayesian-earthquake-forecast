# ------------------------------------------------------------------------------
# 01_model_order_selection.R
#
# This script selects the optimal order `p` for an autoregressive (AR) model 
# applied to annual global earthquake count data (magnitude > 7). It evaluates 
# Partial Autocorrelation Function (PACF), Akaike Information Criterion (AIC), 
# and Bayesian Information Criterion (BIC) for AR orders from 1 to 15.
#
# Output: AIC and BIC plots to help determine the most appropriate AR order.
# ------------------------------------------------------------------------------



earthquakes.dat <- read.delim("data/earthquakes.txt")
earthquakes.dat$Quakes=as.numeric(earthquakes.dat$Quakes)
y=earthquakes.dat$Quakes
##########################################################################
##########################################################################
###################### Determining order of AR by PACF AND  BIC AND DIC PLOT ##########################
#########################################################################
n.all = length(y)
acf(y)
pacf(y, lag.max = 4)


earthquakes.dat <- read.delim("data/earthquakes.txt")
earthquakes.dat$Quakes=as.numeric(earthquakes.dat$Quakes)


y.dat=earthquakes.dat$Quakes#[1:100] ## this is the training data
#y.new=earthquakes.dat$Quakes[101:103] ## this is the test data




n.all=length(y.dat)
p.star=15
Y=matrix(y.dat[(p.star+1):n.all],ncol=1)
sample.all=matrix(y.dat,ncol=1)
n=length(Y)
p=seq(1,p.star,by=1)


design.mtx=function(p_cur){
  Fmtx=matrix(0,ncol=n,nrow=p_cur)
  for (i in 1:p_cur) {
    start.y=p.star+1-i
    end.y=start.y+n-1
    Fmtx[i,]=sample.all[start.y:end.y,1]
  }
  return(Fmtx)
}




criteria.ar=function(p_cur){
  Fmtx=design.mtx(p_cur)
  beta.hat=chol2inv(chol(Fmtx%*%t(Fmtx)))%*%Fmtx%*%Y
  R=t(Y-t(Fmtx)%*%beta.hat)%*%(Y-t(Fmtx)%*%beta.hat)
  sp.square=R/(n-p_cur)
  aic=2*p_cur+n*log(sp.square)
  bic=log(n)*p_cur+n*log(sp.square)
  result=c(aic,bic)
  return(result)
}


criteria=sapply(p,criteria.ar)


plot(p,criteria[1,],type='p',pch='a',col='red',xlab='AR order p',ylab='Criterion',main='',
     ylim=c(min(criteria)-10,max(criteria)+10))
points(p,criteria[2,],pch='b',col='blue')
