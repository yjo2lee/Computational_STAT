# Riemann method
myRiemann<- function(f, a, b, threshold = 10^(-7)){
  h = b-a; Rn = 0
  maxiter = 1000; iter = 0; Rn.0 = 0; error = 1
  while(error>=threshold && iter < maxiter){
    h = h/(2)
    i = seq(0, 2^(iter+1)-1, 1)
    Rn = h*sum(f(a+i*h))
    error = abs(Rn - Rn.0)
    Rn.0 = Rn
    iter = iter + 1
  } 
  return(c(Rn, iter))
}
round(myRiemann(g1, 0,2),4)
g1<- function(x){
  return(x^2)
}
# Trapezoidal rule
myTPz<- function(f, a, b, threshold = 10^(-7)){
  h = b-a; Tn = 0
  maxiter = 100; iter = 0; Tn.0 = 0; error = 1
  while(error>=threshold && iter < maxiter){
    h = h/(2)
    i = seq(1, 2^(iter+1)-1, 1)
    Tn = h*(f(a)/2+sum(f(a+i*h))+f(b)/2)
    error = abs(Tn - Tn.0)
    Tn.0 = Tn
    iter = iter + 1
  } 
  return(c(Tn, iter))
}
round(myTPz(f1, 0,2),4)

mySimpson<- function(f, a, b, threshold = 10^(-7)){
  h = b-a; Sn = 0
  maxiter = 100; iter = 0; Sn.0 = 0; error = 1
  while(error>=threshold && iter < maxiter){
    h = h/(2)
    i = seq(1, 2^(iter), 1)
    Sn = (h/3)*(sum(f(a+(2*i-2)*h)+4*f(a+(2*i-1)*h)+f(a+(2*i)*h)))
    error = abs(Sn - Sn.0)
    Sn.0 = Sn
    iter = iter + 1
  } 
  return(c(Sn, iter))
}
round(mySimpson(f1, 0,2),4)


posterior.1<- function(mu){
  x<-c(6.52, 8.32, 0.31, 2.82, 9.96, 0.14, 9.64)
  return(dcauchy(mu, location = 5, scale = 2) * dnorm(mean(x), mean = mu, sd = 3/sqrt(7)))
}
posterior.2<- function(mu){
  x<-c(6.52, 8.32, 0.31, 2.82, 9.96, 0.14, 9.64)
  return(7.84654*dcauchy(mu, location = 5, scale = 2) * dnorm(mean(x), mean = mu, sd = 3/sqrt(7)))
}
integrate(posterior.1, -Inf, Inf)
round(myRiemann(posterior.1, -10^6, 10^6),4)
round(myTPz(posterior.1, -10^6, 10^6),4)
round(mySimpson(posterior.1, -10^6, 10^6),4)
system.time(myRiemann(posterior.1, -10^6, 10^6))
system.time(myTPz(posterior.1, -10^6, 10^6))
system.time(mySimpson(posterior.1, -10^6, 10^6))
system.time(for(i in 1:10^3) integrate(posterior.1, -Inf, Inf))

integrate(posterior.2, 2, 8)
myRiemann(posterior.2, 2, 8)
myTPz(posterior.2, 2, 8)
mySimpson(posterior.2, 2, 8)
system.time(myRiemann(posterior.2, 2, 8))
system.time(for(i in 1:10^3) myTPz(posterior.2, 2, 8))
system.time(for(i in 1:10^3) mySimpson(posterior.2, 2, 8))
system.time(for(i in 1:10^3) integrate(posterior.2, 2, 8))

u = exp(mu)/(1+exp(mu))

posterior.3<- function(u){
  x<-c(6.52, 8.32, 0.31, 2.82, 9.96, 0.14, 9.64)
  jac = 1/(u*(1-u))
  return(jac*7.84654*dcauchy(log(u/(1-u)), location = 5, scale = 2) * dnorm(mean(x), mean = log(u/(1-u)), sd = 3/sqrt(7)))
}

myRiemann(posterior.3, exp(3)/(1+exp(3)), 1)
myTPz(posterior.3, exp(3)/(1+exp(3)), 0.99999)
mySimpson(posterior.3, exp(3)/(1+exp(3)), 0.99999)
integrate(posterior.3, exp(3)/(1+exp(3)), 1)
system.time(myRiemann(posterior.3, exp(3)/(1+exp(3)), 1))
system.time(for(i in 1:10^3) integrate(posterior.3, exp(3)/(1+exp(3)), 1))

myRiemann(posterior.3, exp(3)/(1+exp(3)), 0.99999)
myTPz(posterior.3, exp(3)/(1+exp(3)), 0.99999)
mySimpson(posterior.3, exp(3)/(1+exp(3)), 0.99999)
integrate(posterior.3, exp(3)/(1+exp(3)), 0.99999)
system.time(myRiemann(posterior.3, exp(3)/(1+exp(3)), 0.99999))
system.time(for(i in 1:10^3) myTPz(posterior.3, exp(3)/(1+exp(3)), 0.99999))
system.time(for(i in 1:10^3) mySimpson(posterior.3, exp(3)/(1+exp(3)), 0.99999))
system.time(for(i in 1:10^3) integrate(posterior.3, exp(3)/(1+exp(3)), 0.99999))

posterior.4<- function(u){
  x<-c(6.52, 8.32, 0.31, 2.82, 9.96, 0.14, 9.64)
  jac = -1/u^2
  return(jac*7.84654*dcauchy(1/u, location = 5, scale = 2) * dnorm(mean(x), mean = 1/u, sd = 3/sqrt(7)))
}
integrate(posterior.4, 1/3, 0)
myRiemann(posterior.4, 1/3, 0)
myTPz(posterior.4, 1/3, 0)
mySimpson(posterior.4, 1/3, 0)
system.time(myRiemann(posterior.4, 1/3, 0))
system.time(for(i in 1:10^3) integrate(posterior.4, 1/3, 0))

integrate(posterior.4, 1/3, 0.000001)
myRiemann(posterior.4, 1/3, 0.000001)
myTPz(posterior.4, 1/3, 0.000001)
mySimpson(posterior.4, 1/3, 0.000001)
system.time(myRiemann(posterior.4, 1/3, 0.000001))
system.time(for(i in 1:10^3) myTPz(posterior.4, 1/3, 0.000001))
system.time(for(i in 1:10^3) mySimpson(posterior.4, 1/3, 0.000001))
system.time(for(i in 1:10^3) integrate(posterior.4, 1/3, 0.000001))

# Romberg's algo
f<- function(x){
  1/x
}

myTn<- function(f=f, a, b, n){
  h = (b-a)/(n); Tn = 0
  if(n == 0)
  i = seq(1, n, 1)
  Tn = h*(f(a)/2+sum(f(a+i*h))+f(b)/2)
  return(Tn)
}
trapij<- function(f=f, a, b, i, j){
  if (j == 0){
    if (i == 0){
      return(2.4)
    }
    else{
      return(myTn(f, a, b, 2^i))
    }
  }
  else{
    return({4^j * trapij(f, 1,5,i,j-1) - trapij(f, 1,5,i-1, j-1)}/(4^j -1))
  }
}

Qij<- function(f,a, b, i,j){
  (trapij(f, a,b,i,j)-trapij(f,a,b,i-1, j))/(trapij(f,a,b,i+1, j)-trapij(f,a,b,i,j))
}

## Rejection Sampling: Gamma Deviates
q<- function(y,a){
  exp(a*log(t(y,a)/a)-t(y,a)+a)
}
e<- function(y){
  exp(-y^2/2)
}
q.e<- function(z,a){
  exp(z^2/2 + a*log(t(z,a)/a)-t(z,a)+a)
}
t<- function(z, a=a){
  b = 1/sqrt(9*a)
  a*(1+b*z)^3
}

GammaDer<- function(r){
  n = 100000; a = r-1/3; b = 1/sqrt(9*a)
  z = rnorm(n,0,1) ; z = z[z>-1/b] ; u = runif(length(z),0,1)
  raw.x = t(z, a)
  accepted = raw.x[t(z,a) > 0 & u <= q.e(z,a)]
  rate = length(accepted)/length(z)
  par(mfrow=c(1,2))
  x = seq(0, 10, 0.01)
  plot(x, e(x), type="l", main = paste("target and envelope: r =", r))
  lines(x, q(x, a), type = 'l', col = 2)
  legend(4,1,c("Envelope","Target"),lty=2:1,lwd=2:1,col=1:2) 
  qqplot(qgamma((1:n)/(n+1),r,1), accepted, xlab="Theoretical Quantiles", ylab="Sample Quantiles",
         main=paste("Q-Q Plot: r=",r)) ; abline(a=0,b=1,col=2)
  par(mfrow=c(1,1))
  return(c(rate, 1-rate, mean(accepted), sd(accepted)))
}

GammaDer(1); GammaDer(2); GammaDer(3); GammaDer(4); GammaDer(5)

# Rejection Sampling a Bayesian Posterior
q.lam<- function(l){
  x<- c(8, 3, 4, 3, 1, 7, 2, 6, 2, 7)
  return(prod(dpois(x, l)))
}
x<- c(8, 3, 4, 3, 1, 7, 2, 6, 2, 7)
Bay.posterior<- function(x){
  n = 10000;
  l = rlnorm(n, log(4), 0.5) ; u = runif(n,0,1)
  q = sapply(l,q.lam)
  e = prod(dpois(x, mean(x)))
  accepted = l[u < q/e]
  rate = length(accepted)/n
  y = seq(0, 20, 0.01)
  plot(y, dlnorm(y, log(4), 0.5)*prod(dpois(x, mean(x))), type="l", main = paste("target and envelope"), ylab = 'density')
  lines(y, dlnorm(y, log(4), 0.5)*sapply(y,q.lam), type = 'l', col = 2)
  legend(4,1,c("Envelope","Target"),lty=2:1,lwd=2:1,col=1:2) 
  return(c(rate, 1-rate, mean(accepted), sd(accepted)))
}

Bay.posterior(x)






 