require("numDeriv")

g1<-function(x){
  log(x)/(1+x)
}

f3<-function(x){
  (1+1/x-log(x))/(1+x)^2
}

g2 <- function(x){
   3^x/sqrt(x)
}

g3 <- function(x){
  cos(x+1/x)/sqrt(x)
}

# Define function
myBisec<-function(f, a, b, threshold = 10^-8){
  maxiter<-1000
  err<-1
  niter<-0
  x0<-(a+b)/2
  while ( niter<=maxiter && err >= threshold){
    #update interval
    if (grad(f, a)*grad(f, x0) <=0)	{b<-x0}
    else	{a<-x0}
    #update x
    oldx0<-x0
    x0<-(a+b)/2
    #update error and niter
    err<-abs(oldx0-x0)
    niter<-niter+1
    #print(paste(niter,x0,err,sep="  "))
  }
  df = data.frame(algorithm = "bisection", initial = c(a,b), x_star = x0, optimag = f(x0), niter = niter)
  return(df)
}

myNewton<-function(f, x0, threshold = 10^-8){
  maxiter = 1000
  err = 1
  niter = 0
  x.init = x0
  while (niter<=maxiter && err >= threshold){
    #update x
    x1 = x0 - (grad(f, x0)/hessian(f, x0)[1])
    err = abs(x0-x1)
    x0 = x1
    niter = niter+1
    #print(paste(niter,x0,err,sep="  "))
  }
  df = data.frame(algorithm = "Newton", initial = x.init, x_star = round(x1, 4), optimag = round(f(x1), 4), niter = niter)
  return(df)
}

mySecant<-function(f, x0, x1, threshold = 10^-8){
  maxiter = 1000
  err = 1
  niter = 0
  x0.init = x0; x1.init = x1
  while (niter<=maxiter && err >= threshold){
    #update x
    x2 = x1 - grad(f, x1)*(x1-x0)/(grad(f,x1)-grad(f,x0))
    #update error and niter
    err = abs(x2-x1)
    x0 = x1
    x1 = x2
    niter = niter+1
    # to see each step
    #print(paste(niter,x1,err,sep="  "))
  }
  df = data.frame(algorithm = "Secant", initial = c(x0.init, x1.init), x_star = round(x2,4), optimag = round(f(x2),4), niter = niter)
  return(df)
}

# result
plot(seq(1,20,0.1),g1(seq(1,20,0.1)),type='l',xlab="x",ylab="g1(x)",main="g1")
myBisec(g1, 1,5)
myNewton(g1, 1)
mySecant(g1, 1, 5)
optimize(g1,c(1,5), maximum = T)
system.time(for(i in 1:10^3) myBisec(g1, 1,5))
system.time(for(i in 1:10^3) myNewton(g1, 1))
system.time(for(i in 1:10^3) mySecant(g1, 1,5))
system.time(for(i in 1:10^3) optimize(g1,c(1,5), maximum = T))


plot(seq(1,10,0.1),g2(seq(1,10,0.1)),type='l',xlab="x",ylab="g2(x)",main="g2")
myBisec(g2, 1,4)
myNewton(g2, 1)
mySecant(g2, 1, 4)
system.time(for(i in 1:10^3) myBisec(g2, 1,5))
system.time(for(i in 1:10^3) myNewton(g2, 1))
system.time(for(i in 1:10^3) mySecant(g2, 1,5))
system.time(for(i in 1:10^3) optimize(g2,c(1,5), maximum = T))

plot(seq(1,10,0.1),g3(seq(1,10,0.1)),type='l',xlab="x",ylab="g3(x)",main="g3")
myBisec(g3, 5,7)
myNewton(g3, 5)
mySecant(g3, 5, 7)
system.time(for(i in 1:10^3) myBisec(g3, 5,7))
system.time(for(i in 1:10^3) myNewton(g3, 5))
system.time(for(i in 1:10^3) mySecant(g3, 5,7))
system.time(for(i in 1:10^3) optimize(g3,c(5,7), maximum = T))


#Exercise2.1
#a) log-likelihood
x.prob1 = c(1.77,-0.23,2.76,3.80,3.47,56.75,-1.34,4.24,-2.44,3.29,3.71,-2.40,
  4.53,-0.07,-1.05,-13.87,-2.53,-1.75,0.27,43.21)
loglkd.cauchy <- function(theta){
  loglkd.theta = rep(0, length(theta))
  for(i in 1:length(theta)){
    loglkd.theta[i] = sum(dcauchy(x.prob1, location = theta[i], scale = 1, log = T))
  }
  return(loglkd.theta)
}
plot(seq(-15,15,0.5), loglkd.cauchy(seq(-15,15,0.5)), type = 'l', xlab = 'theta', ylab='loglikelihood', main = "loglikelihood of cauchy(theta, 1)")
myNewton(loglkd.cauchy, 0)
optimize(loglkd.cauchy, c(-3,1), maximum = T)
theta.1<-c(-11,-1,0,1.5,4,4.7,7,8,38)

for(i in 2:length(theta.1)){
  print(myNewton(loglkd.cauchy, theta.1[i]))
}
#b)
myBisec(loglkd.cauchy, -1, 1)
myBisec(loglkd.cauchy, -3, -1)
myBisec(loglkd.cauchy, -2, 0)
myBisec(loglkd.cauchy, 0, 2)
myBisec(loglkd.cauchy, 1, 3)
myBisec(loglkd.cauchy, -2, 2)
myBisec(loglkd.cauchy, -3, 3)

#d)
mySecant(loglkd.cauchy, -2, -1)
mySecant(loglkd.cauchy, -3, 3)
mySecant(loglkd.cauchy, -1, 1)
mySecant(loglkd.cauchy, 0, 2)
mySecant(loglkd.cauchy, 1, 3)
mySecant(loglkd.cauchy, -2, 2)
mySecant(loglkd.cauchy, -2, 0)

#e)
system.time(for(i in 1:10^3) myBisec(loglkd.cauchy, -1,1))
system.time(for(i in 1:10^3) myNewton(loglkd.cauchy, 0))
system.time(for(i in 1:10^3) mySecant(loglkd.cauchy, -2,1))
system.time(for(i in 1:10^3) optimize(loglkd.cauchy, c(-3,1), maximum = T))

#a_normal)
x.prob.norm = rnorm(20, 0, 1)
loglkd.normal <- function(theta){
  loglkd.theta = rep(0, length(theta))
  for(i in 1:length(theta)){
    loglkd.theta[i] = sum(dnorm(x.prob.norm, mean = theta[i], sd = 1, log = T))
  }
  return(loglkd.theta)
}
for(i in 2:length(theta.1)){
  print(myNewton(loglkd.normal, theta.1[i]))
}

#b_normal)
myBisec(loglkd.normal, -1, 1)
myBisec(loglkd.normal, -3, -1)
myBisec(loglkd.normal, -2, 0)
myBisec(loglkd.normal, -2, 2)
myBisec(loglkd.normal, 1, 3)


#d_normal)
mySecant(loglkd.normal, -1, -1)
mySecant(loglkd.normal, -3, -1)
mySecant(loglkd.normal, 0, 2)
mySecant(loglkd.normal, -2, 2)
mySecant(loglkd.normal, -2, 0)

optimize(loglkd.normal, c(-1,1), maximum = T)

#e)
system.time(for(i in 1:10^3) myBisec(loglkd.normal, -1,1))
system.time(for(i in 1:10^3) myNewton(loglkd.normal, -1))
system.time(for(i in 1:10^3) mySecant(loglkd.normal, -2,0))
system.time(for(i in 1:10^3) optimize(loglkd.normal, c(-1,1), maximum = T))

# 2.2

# optim
# default: Nelder-Mead"
f <- expression(x^2, 'x')

# symbolic diffrentiation
D(f, 'x')
D(D(f, 'x'),'x')

# symbolic integration; maple