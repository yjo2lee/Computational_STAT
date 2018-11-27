library(gplots)

# Option pricing
## arithmetic mean
mc.ari.array = rep(0,n)
sd.array = rep(0,100)
for(j in 1:100){
  Asian.opt.avg = 0
  for (i in 1:n){
    t = 50
    z.t = rnorm(t, 0, 1)
    r = 0.05; vol = 0.3; s.init = 100; k = 102
    s.t = s.init*cumprod(exp((r-vol^2/2)/365 + vol*z.t/sqrt(365)))
    Asian.opt = exp(-r*t/365)*max(0, mean(s.t)-k)
    Asian.opt.avg = Asian.opt.avg + Asian.opt/n
    mc.ari.array[i] = Asian.opt
    i = i+1
  }
  sd.array[j] = Asian.opt.avg
}

for (i in 1:n){
  t = 50
  z.t = rnorm(t, 0, 1)
  r = 0.05; vol = 0.3; s.init = 100; k = 102
  s.t = s.init*cumprod(exp((r-vol^2/2)/365 + vol*z.t/sqrt(365)))
  Asian.opt = exp(-r*t/365)*max(0, mean(s.t)-k)
  Asian.opt.avg = Asian.opt.avg + Asian.opt/n
  mc.ari.array[i] = Asian.opt
  i = i+1
}
mc.ari.sd = sd(mc.ari.array)/sqrt(n)

## geometric_mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
Asian.opt.avg.geo = 0
mc.geo.array = rep(0,n)
for (i in 1:n){
  t = 50
  z.t = rnorm(t, 0, 1)
  r = 0.05; vol = 0.3; s.init = 100; k = 102
  s.t = s.init*cumprod(exp((r-vol^2/2)/365 + vol*z.t/sqrt(365)))
  Asian.opt = exp(-r*t/365)*max(0, gm_mean(s.t)-k)
  Asian.opt.avg.geo = Asian.opt.avg.geo + Asian.opt/n
  mc.geo.array[i] = Asian.opt
  i = i+1
}
geo_array = rep(0,100)
for(j in 1:100){
  Asian.opt.avg.geo = 0
  for (i in 1:n){
    t = 50
    z.t = rnorm(t, 0, 1)
    r = 0.05; vol = 0.3; s.init = 100; k = 102
    s.t = s.init*cumprod(exp((r-vol^2/2)/365 + vol*z.t/sqrt(365)))
    Asian.opt = exp(-r*t/365)*max(0, gm_mean(s.t)-k)
    Asian.opt.avg.geo = Asian.opt.avg.geo + Asian.opt/n
    i = i+1
  }
  geo_array[j] = Asian.opt.avg.geo
}

mc.geo.sd = sd(mc.geo.array)/sqrt(n)
mu.cv = function(mu.mc, lambda, theta.mc,theta){
  return(mu.mc + lambda*(theta.mc-theta))
}

## control_variate method
mu.cv.array = rep(0, 100)
for(i in 1:100){
  c3 = 1+1/t
  c2 = vol*sqrt(c3*t/1095*(1+1/2/t))
  c1 = 1/c2*(log(s.init/k)+c3*t/730*(r-vol^2/2)+c3*vol^2*t/1095*(1+1/2/t))
  theta.cv = s.init*pnorm(c1)*exp(-t*(r+c3*vol^2/6)*(1-1/t)/730)-k*pnorm(c1-c2)*exp(-r*t/365)
  asian.mu.cv = mu.cv(sd.array[i], -1.0217, geo_array[i], theta.cv)
  mu.cv.array[i] = asian.mu.cv
  i = i+1
}
sd.cv = sd(mu.cv.array)/sqrt(100)
lambda.hat = -cov(sd.array, geo_array)/var(geo_array)

# 6.3
#importance sampling with std
##target
q<-function(x){ 
  return(exp(-abs(x)^3/3)) 
}
##envelope
e<-function(x){ 
  return(dnorm(x,0,1))
}
##h(x)
h<-function(x) {  
  return(x^2)
}
x<-rnorm(n,0,1) 
w<-q(x)/e(x); w<-w/sum(w)
mu.std.is = sum(h(x)*w)

curve(e(x),-3:3,xlim=c(-3,3),ylim=c(0,1),lwd=2,xlab="x",ylab="",col="darkblue") 
curve(q(x),-3:3,xlim=c(-3,3),col=2,lty=2,lwd=2,add=TRUE)
legend('topright' , c('envelope' , 'target') , col=c('black','red'), lty=1:2,lwd=2)

#rejection sampling
u<-runif(n,0,1)
y<-x[u<=q(x)/e(x)/sqrt(2*pi)] 
c(Accept=length(y)/10^5,Reject=1-length(y)/10^5) 
mean(h(y))
curve(e(x)*sqrt(2*pi),-3:3,xlim=c(-3,3),ylim=c(0,1),lwd=2,xlab="x",ylab="",col="darkblue") 
curve(q(x),-3:3,xlim=c(-3,3),col=2,lty=2,lwd=2,add=TRUE)
legend('topright' , c('envelope' , 'target') , col=c('darkblue','red'), lty=1:2,lwd=2)

# importance-weight: employ a Riemann sum strategy with random nodes
y.order<-sort(y) 
num<-rep(0,length(y)-1) 
denom<-num
for(i in 1:(length(y)-1)){
  num[i]<-(y.order[i+1]-y.order[i])*h(y.order[i])*q(y.order[i])
  denom[i]<-(y.order[i+1]-y.order[i])*q(y.order[i]) }
sum(num)/sum(denom) #sigma2

# compare the variances
result<-matrix(nrow=100,ncol=2) 
colnames(result)<-c("b","c")
for(j in 1:100){ 
  x<-rnorm(10^4,0,1)
#Rejection Sampling 
  u<-runif(10^4,0,1) 
  y<-x[u<=q(x)/e(x)/sqrt(2*pi)] 
  result[j,1]<-mean(h(y))
#Riemann summ strategy 
  y.order<-sort(y) 
  num<-rep(0,length(y)-1) 
  denom<-num
  for(i in 1:(length(y)-1)){
    num[i]<-(y.order[i+1]-y.order[i])*h(y.order[i])*q(y.order[i])
    denom[i]<-(y.order[i+1]-y.order[i])*q(y.order[i]) 
  }
    result[j,2]<-sum(num)/sum(denom) 
}
mu<-apply(result,2,mean) 
Var<-apply(result,2,var)

# 6.6: Testing the hypotheses H0: lambda = 2
## alpha: type-I-error: p(reject H0|H0 is true)
z = function(x){
  (mean(x)-2)/sqrt(2/25)
}

h <- function(x){ifelse(x > 1.645, 1, 0)}
n2 = 5000

##naive mc
mc.alpha = rep(0, n2)
for (i in 1:n2){
  x = rpois(25, 2)
  mc.alpha[i] = mean(sapply(z(x), h))
}
x = rpois(25, 2); mean(sapply(z(x), h))
mean(mc.alpha); sd(mc.alpha)/sqrt(n2)

##antithetic
as.alpha = rep(0, n2)
for (i in 1:n2){
  u1 = runif(25, 0,1); u2 = 1-u1
  x1 = qpois(u1, 2); x2 = qpois(u2, 2)
  as.alpha[i] = (mean(sapply(z(x1), h)) + mean(sapply(z(x2), h)))/2
}
mean(as.alpha); sd(as.alpha)/sqrt(n2)


## importance sampling with unstd weights & std weights
q.p = function(x){dpois(x, 2*25)}
e.p = function(x){dpois(x,2.4653*25)}

h2 = function(x){ifelse(x >= 1.645*sqrt(25*2)+25*2, 1, 0)}
unstd.mean = rep(0, 100); std.mean = rep(0, 100)
for(j in 1:100){
  is.unstd = rep(0, n2)
  weightsum = 0
  for (i in 1:n2){
    x = rpois(25, 2.4653) 
    y = sum(x)
    w = q.p(y)/e.p(y)
    is.unstd[i] = h2(y)*w
    weightsum = weightsum+w
    #is.std[i] = h2(y)*w.s
  }
  unstd.mean[j] = mean(is.unstd)
  std.mean[j] = sum(is.unstd/weightsum)
}

## control variate for importance sampling
cv.is.mean = rep(0, 100)
for (j in 1:100){
  cv.is = rep(0, n2); weight = rep(0, n2); h3 = rep(0, n2)
  for (i in 1:n2){
    x = rpois(25, 2.4653) 
    y = sum(x)
    w = q.p(y)/e.p(y)
    weight[i] = w
    h3[i] = h2(y)
  }
  lambda = -lm(h3*weight ~ weight)$coef[2]
  cv.is = mean(is.unstd) + lambda*(mean(weight)-1)
  cv.is.mean[j] = mean(cv.is)
}

makeCI = function(mu, sd){
  return(c(round(mu-1.96*sd, 4), round(mu + 1.96*sd,4)))
}

## power: p(accept H1|H1 is true)
lam = seq(2.2, 4, length = 5)
z = function(x){
  (mean(x)-2)/sqrt(2/25)
}
hb <- function(x){ifelse(x < 1.645, 1, 0)}
n2 = 5000

##naive mc
naive_power = function(lam){
  mc.power = rep(0, n2)
  for (i in 1:n2){
    x = rpois(25, lam)
    mc.power[i] = mean(sapply(z(x), h))
  }
  return(c(mean(mc.power), sd(mc.power)/sqrt(n2)))
}
mc.res = sapply(lam, naive_power)

## antithetic
as_power = function(lam){
  as.power = rep(0, n2)
  for (i in 1:n2){
    u1 = runif(25,0,1); u2 = 1-u1
    x1 = qpois(u1, lam); x2 = qpois(u2, lam)
    as.power[i] = (mean(sapply(z(x1), h)) + mean(sapply(z(x2), h)))/2
  }
  return(c(mean(as.power), sd(as.power)/sqrt(n2)))
}
as.res = sapply(lam, as_power)

## IS
IS_power = function(lam){
  q.p = function(x){dpois(x, lam*25)}
  e.p = function(x){dpois(x,(lam+0.4653)*25)}
  h2 = function(x){ifelse(x >= 1.645*sqrt(25*2)+25*2, 1, 0)}
  unstd.mean = rep(0, 100); std.mean = rep(0, 100)
  for(j in 1:100){
    is.unstd = rep(0, n2)
    weightsum = 0
    for (i in 1:n2){
      x = rpois(25, (lam + 0.4653)) 
      y = sum(x)
      w = q.p(y)/e.p(y)
      is.unstd[i] = h2(y)*w
      weightsum = weightsum+w
      #is.std[i] = h2(y)*w.s
    }
    unstd.mean[j] = mean(is.unstd)
    std.mean[j] = sum(is.unstd/weightsum)
  }
  return(c(mean(unstd.mean), sd(unstd.mean), mean(std.mean), sd(std.mean)))
}
IS_res = sapply(lam, IS_power)

## control variate
h2 = function(x){ifelse(x >= 1.645*sqrt(25*2)+25*2, 1, 0)}
cv_power = function(lam){
  q.p = function(x){dpois(x, lam*25)}
  e.p = function(x){dpois(x,(lam+0.4653)*25)}
  cv.is.mean = rep(0, 100)
  for (j in 1:100){
    cv.is = rep(0, n2); weight = rep(0, n2); h3 = rep(0, n2)
    for (i in 1:n2){
      x = rpois(25, (lam+0.4653)) 
      y = sum(x)
      w = q.p(y)/e.p(y)
      weight[i] = w
      h3[i] = h2(y)
    }
    lambda = -lm(h3*weight ~ weight)$coef[2]
    cv.is = mean(h3*weight) + lambda*(mean(weight)-1)
    cv.is.mean[j] = cv.is
  }
  return(c(mean(cv.is.mean), sd(cv.is.mean)))
}
cv_res = sapply(lam, cv_power)

library(gplots) 
plotCI(lam,mc.res[1,],uiw=1.96*mc.res[2,],liw=1.96*mc.res[2,],pch=20,
         gap=0,type="o",cex=0.95,xlab=expression(lambda),ylab="power",main="Power Curve of Naive MC")
plotCI(lam,as.res[1,],uiw=1.96*as.res[2,],liw=1.96*as.res[2,],pch=20,
       gap=0,type="o",cex=0.95,xlab=expression(lambda),ylab="power",main="Power Curve of Antithetic sampling")
plotCI(lam,IS_res[1,],uiw=1.96*IS_res[2,],liw=1.96*IS_res[2,],pch=20,
       gap=0,type="o",cex=0.95,xlab=expression(lambda),ylab="power",main="Power Curve of IS(unstd.)")
plotCI(lam,IS_res[3,],uiw=1.96*IS_res[4,],liw=1.96*IS_res[4,],pch=20,
       gap=0,type="o",cex=0.95,xlab=expression(lambda),ylab="power",main="Power Curve of IS(std.)")
plotCI(lam,cv_res[1,],uiw=1.96*cv_res[2,],liw=1.96*cv_res[2,],pch=20,
       gap=0,type="o",cex=0.95,xlab=expression(lambda),ylab="power",main="Power Curve of control variate")

write.csv(result,"EX_6_6_b_result.csv")


# 6.7 pricing a European call option
## mc
n3 = 10^5
z.t = rnorm(n3,0,1)
s.t = s.init*exp((r-vol^2/2)*t/365+vol*z.t*sqrt(t/365)) 
c1 = exp(-r*t/365)*sapply(s.t-k,function(x){max(0,x)}) 
mc.mean.a = mean(c1)
mc.sd.a = sd(c1)/sqrt(n3)

## Balck-Scholes
t.s<-t/365; w<-(r*t.s+ vol^2*t.s/2-log(k/s.init))/vol/sqrt(t.s) 
c2<-s.init*pnorm(w)-k*exp(-r*t.s)*pnorm(w-vol*sqrt(t.s)) 
c2

## asian call option: mc
mc.opt.array = rep(0,100)
for(j in 1:100){
  opt.avg = 0
  for (i in 1:n3){
    t = 30
    z.t = rnorm(t, 0, 1)
    r = 0.05; vol = 0.5; s.init = 50; k = 52
    s.t = s.init*cumprod(exp((r-vol^2/2)/365 + vol*z.t/sqrt(365)))
    opt = exp(-r*t/365)*max(0, mean(s.t)-k)
    opt.avg = opt.avg + opt/n3
    i = i+1
  }
  mc.opt.array[j] = opt.avg
}

## asian call option: cv
geo.opt = rep(0,100)
for(j in 1:100){
  opt.geo.avg = 0
  for (i in 1:n3){
    t = 30
    z.t = rnorm(t, 0, 1)
    r = 0.05; vol = 0.5; s.init = 50; k = 52
    s.t = s.init*cumprod(exp((r-vol^2/2)/365 + vol*z.t/sqrt(365)))
    opt.geo = exp(-r*t/365)*max(0, gm_mean(s.t)-k)
    opt.geo.avg = opt.geo.avg + opt.geo/n3
    i = i+1
  }
  geo.opt[j] = opt.geo.avg
}
mean(geo.opt)
geo.sd = sd(geo.opt)/sqrt(n3)
mu.cv = function(mu.mc, lambda, theta.mc,theta){
  return(mu.mc + lambda*(theta.mc-theta))
}

cv.opt = rep(0, 100)
for(i in 1:100){
  c3 = 1+1/t
  c2 = vol*sqrt(c3*t/1095*(1+1/2/t))
  c1 = 1/c2*(log(s.init/k)+c3*t/730*(r-vol^2/2)+c3*vol^2*t/1095*(1+1/2/t))
  theta.cv = s.init*pnorm(c1)*exp(-t*(r+c3*vol^2/6)*(1-1/t)/730)-k*pnorm(c1-c2)*exp(-r*t/365)
  asian.mu.cv = mu.cv(mc.opt.array[i], -1, geo.opt[i], theta.cv)
  cv.opt[i] = asian.mu.cv
  i = i+1
}
sd.cv = sd(cv.opt)/sqrt(100)
lambda.hat.c = -cov(mc.opt.array, geo.opt)/var(geo.opt)

## antithetic
for(i in 1:(5*10^4)){
  z.t = rnorm(t,0,1) 
  s.as = matrix(NA,ncol=2,nrow=t+1) 
  a.as<-matrix(NA,nrow=5*10^4,ncol=2)
  s.as[1,] = s.init
  for(j in 1:t){
    s.as[j+1,1] = s.as[j,1]*exp((r-vol^2/2)/365+vol*z.t[j]/sqrt(365))
    s.as[j+1,2] = s.as[j,2]*exp((r-vol^2/2)/365+vol*(-z.t[j])/sqrt(365)) 
  }
  a.as[i,1] = exp(-r*t/365)*max(0,mean(s.as[-1,1])-k) 
  a.as[i,2] = exp(-r*t/365)*max(0,mean(s.as[-1,2])-k)
}
