# Part1
# GMM: estimate mixing proportion using independent chain 
d = rbinom(100,1,0.7)
y = d*rnorm(100,7,0.5)+(1-d)*rnorm(100,10,0.5) 
hist(y,prob=T,nclass = 20,col="gray")
rug(y) ; lines(density(y),lwd=2)

#MHRatio
f = function(d){prod(d*dnorm(y,7, 0.5) + (1-d)*dnorm(y, 10,0.5))}  # P(d|y)
g = function(d){dbeta(d, 1, 1)}
R = function(d, d.1){f(d.1)*g(d)/f(d)/g(d.1)}

#Beta(1,1)=Unif(0,1)
n = 10^4
delta1=rep(NA,n) ; delta1[1]=rbeta(1,1,1) 
for(i in 1:(n-1)){
  d=rbeta(1,1,1)
  MH.R=R(delta1[i],d); prob=min(MH.R,1) 
  delta1[i+1]=ifelse(rbinom(1,1,prob)==1,d,delta1[i])
}
plot(delta1,type="l",ylim=c(0,1),ylab="delta",xlab="t",main="A Sample Path of delta") 
summary(delta1) ; sd(delta1)

#Beta(2,10)
delta2=rep(NA,n) ; delta2[1]=rbeta(1,2,10) 
for(i in 1:(n-1)){
  d=rbeta(1,2,10)
  MH.R=R(delta2[i],d)
  prob=min(MH.R,1) 
  delta2[i+1]=ifelse(rbinom(1,1,prob)==1,d,delta2[i])
}
plot(delta2,type="l",ylim=c(0,1),ylab="delta",xlab="t",main="A Sample Path of delta") 
summary(delta2) ; sd(delta2)

# Part2
# GMM: estimate mixing proportion using random walk chain
logit.inverse = function(x){exp(x)/(1+exp(x))}
jac = function(x){(1+exp(x))^2/exp(x)}

f = function(x){prod(logit.inverse(x)*dnorm(y,7, 0.5) + (1-logit.inverse(x))*dnorm(y, 10,0.5))}  # P(d|y)
g = function(x, b){dunif(x, -b, b)}
R = function(x.t, x.star,b){f(x.star)*jac(x.star)/f(x.t)/jac(x.t)}

delta3 = rep(NA,n); delta3[1] = runif(1,0,1); b = 1; u = rep(NA, n); u[1] = runif(1,-1,1)
for(i in 1:(n-1)){
  u.star = u[i] + runif(1, -b, b)
  MH.R = R(u[i], u.star, b); prob=min(MH.R,1)
  u[i+1] = ifelse(rbinom(1,1,prob) == 1, u.star, u[i])
  delta3[i+1] = logit.inverse(u[i+1])
}
plot(delta3,type="l",ylim=c(0,1),ylab="delta",xlab="t",main="A Sample Paths of delta") 
hist(delta3[-(1:200)],breaks=20,col="gray",main="Histogram of delta",xlab="delta") 
summary(delta3) ; sd(delta3)

delta4 = rep(NA,n); delta4[1] = runif(1,0,1); b = 0.01; u = rep(NA, n); u[1] = runif(1,-1,1)
for(i in 1:(n-1)){
  u.star = u[i] + runif(1, -b, b)
  MH.R = R(u[i], u.star, b); prob=min(MH.R,1)
  u[i+1] = ifelse(rbinom(1,1,prob) == 1, u.star, u[i])
  delta4[i+1] = logit.inverse(u[i+1])
}
plot(delta4,type="l",ylim=c(0,1),ylab="delta",xlab="t",main="A Sample Paths of delta") 
hist(delta4[-(1:200)],breaks=20,col="gray",main="Histogram of delta",xlab=expression(delta)) 
summary(delta4) ; sd(delta4)


# PART 3; Prob)7.5 
HT=c(2,4,6,9,9,9,13,14,18,23,31,32,33,34,43) 
HT.censor=c(10,14,14,16,17,18,18,19,20,20,21,21,
             23,24,29,29,30,30,31,31,31,33,35,37,
             40,41,42,42,44,46,48,49,51,53,54,54,55,56) 
HT.delta=c(rep(1,length(HT)),rep(0,length(HT.censor)))
con=c(1,4,6,7,13,24,25,35,35,39) 
con.censor=c(1,1,3,4,5,8,10,11,13,14,14,15,
              17,19,20,22,24,24,24,25,26,26,26,28,
              29,29,32,35,38,39,40,41,44,45,47,47,47,50,50,51) 
con.delta=c(rep(1,length(con)),rep(0,length(con.censor)))

time=c(HT,HT.censor,con,con.censor)
delta=c(HT.delta,con.delta)
group=c(rep(1,length(HT.delta)),rep(2,length(con.delta))) #1=Hormone Treated, 2=Control 
data=cbind(time,delta,group) ; data=as.data.frame(data)

#a) 
result.a=matrix(ncol=7,nrow=4) 
rownames(result.a)=c("HT","HT.c","con","con.c") 
colnames(result.a)=c("Min","Q1","Median","Mean","Q3","Max","Sd") 
result.a[1,]=c(summary(HT),sd=sd(HT)) 
result.a[2,]=c(summary(HT.censor),sd=sd(HT.censor)) 
result.a[3,]=c(summary(con),sd=sd(con)) 
result.a[4,]=c(summary(con.censor),sd=sd(con.censor))
par(mfrow=c(2,2))
hist(HT,xlab="Recurrence Times",main="Hormone Treated",col="gray") 
hist(HT.censor,xlab="Censoring Times",main="Hormone Treated(censored)",col="gray") 
hist(con,xlab="Recurrence Times",main="Control",col="gray") 
hist(con.censor,xlab="Censoring Times",main="Control(censored)",col="gray")

# library(survival)
# my.fit=survfit(Surv(time,delta)~group,data=data)
# par(mfrow=c(1,1)) 
# plot(my.fit,lty=1:2,lwd=1:2,col=1:2,xlab="time",main="Survival Function") 
# legend(30,0.4,c("Hormone Treated","Control"),lty=1:2,lwd=1:2,col=1:2)

#### (c) #### 
a=3;b=1;c=60;d=120 

Gibbs_prob5=function(a,b,c,d,data){
  n=10^4
  theta=rep(NA,n+1) ; theta[1]=0
  tau=rep(NA,n+1) ; tau[1]=1
  for(i in 1:n){
    theta[i+1]=rgamma(1,shape=sum(data$delta)+a+1,rate=sum(data$time[data$group==2])+c+tau[i]*sum(data$time[data$group==1])+tau[i]*d)
    tau[i+1]=rgamma(1,shape=sum(data$delta[data$group==1])+b+1,rate=theta[i+1]*sum(data$time[data $group==1])+theta[i+1]*d)
  }
  return(list(theta=theta,tau=tau)) 
}
result.c = Gibbs_prob5(a,b,c,d,data) 
plot(result.c$theta[-(1:500)],type="l",xlab="t",ylab=expression(theta),main="A Sample Path of theta") 
plot(result.c$tau[-(1:500)],type="l",xlab="t",ylab=expression(tau),main="A Sample Path of tau") 
acf(result.c$theta[-(1:500)],main="ACF for theta")
acf(result.c$tau[-(1:500)],main="ACF for tau")

#d) 
post=result.c$theta[-(1:500)]^(sum(data$delta)+a)*result.c$tau[-(1:500)]^(sum(data$delta[data$group==1])+b)*exp(-result.c$theta[-(1:500)]*(sum(data$time[data$group==2])+c)-result.c$tau[-(1:500)]*result.c$theta[-(1:500)]*(sum(data$delta[data$group==1])+d))
summary(post)
result.d=matrix(ncol=4,nrow=2)
rownames(result.d)=c("theta","tau")
colnames(result.d)=c("mean","sd","95%lower","95%upper") 
result.d[1,]=c(mean(result.c$theta[-(1:500)]),sd(result.c$theta[-(1:500)]),quantile(result.c$theta[-(1:500)],c(0.25,0.75))) 
result.d[2,]=c(mean(result.c$tau[-(1:500)]),sd(result.c$tau[-(1:500)]),quantile(result.c$tau[-(1:500)],c(0.25,0.75)))

#e) 
tau.prior=function(tau,a,b,c,d){tau^b*gamma(a+1)/(c+d*tau)^(a+1)}
x=seq(0,6,by=0.01) 
plot(x,tau.prior(x,a,b,c,d)/integrate(tau.prior,0,Inf,a=a,b=b,c=c,d=d)$value,type="l",xlim=c(0,6),ylab="density",xlab=expression(tau),main="density of tau")
lines(density(result.c$tau),lwd=2,lty=2,col=2)
legend(4.5,1,c("Prior","Posterior"),lwd=1:2,lty=1:2,col=1:2)

#f)
t.test(result.c$tau, mu=1,alternative="two.sided",conf.level=0.95) 

#g)
##g)_1
##c)
a=1.5;b=0.5;c=30;d=60
result.g1=Gibbs_prob5(a,b,c,d,data) 
par(mfrow=c(2,2))
plot(result.g1$theta[-(1:500)],type="l",xlab="t",ylab=expression(theta),main="A Sample Path of theta") 
plot(result.g1$tau[-(1:500)],type="l",xlab="t",ylab=expression(tau),main="Sample Path of tau") 
acf(result.g1$theta[-(1:500)],main="ACF of theta")
acf(result.g1$tau[-(1:500)],main="ACF of tau")
##d)
post=result.g1$theta[-(1:500)]^(sum(data$delta)+a)*result.g1$tau[-(1:500)]^(sum(data$delta[data$group==1])+b)*exp(-result.g1$theta[-(1:500)]*(sum(data$time[data$group==2])+c)-result.g1$tau[-(1:500)]*result.g1$theta[-(1:500)]*(sum(data$delta[data$group==1])+d))
summary(post)
result.d1=matrix(ncol=4,nrow=2)
rownames(result.d1)=c("theta","tau")
colnames(result.d1)=c("mean","sd","95%lower","95%upper") 
result.d1[1,]=c(mean(result.g1$theta[-(1:500)]),sd(result.g1$theta[-(1:500)]),quantile(result.g1$theta[-(1:500)],c(0.25,0.75))) 
result.d1[2,]=c(mean(result.g1$tau[-(1:500)]),sd(result.g1$tau[-(1:500)]),quantile(result.g1$tau[-(1:500)],c( 0.25,0.75)))
##e)
par(mfrow=c(1,1))
x=seq(0,6,by=0.01) 
plot(x,tau.prior(x,a,b,c,d)/integrate(tau.prior,0,Inf,a=a,b=b,c=c,d=d)$value,ylab="density",xlab=expression(tau),type="l",xlim=c(0,6),ylim=c(0,1),main="density of tau") 
lines(density(result.g1$tau),lwd=2,lty=2,col=2)
legend(4,1,c("Prior","Posterior"),lwd=1:2,lty=1:2,col=1:2)
##f)
t.test(result.g1$tau,mu=1,alternative="greater",conf.level=0.95)

##g)-2
##c
a=6;b=2;c=60;d=240
result.g2=Gibbs_prob5(a,b,c,d,data) 
par(mfrow=c(2,2))
plot(result.g2$theta[-(1:500)],type="l",xlab="t",ylab=expression(theta),main="A Sample Path of theta") 
plot(result.g2$tau[-(1:500)],type="l",xlab="t",ylab=expression(tau),main="A Sample Path of tau") 
acf(result.g2$theta[-(1:500)],main="ACF of theta")
acf(result.g2$tau[-(1:500)],main="ACF of tau")
##d 
post=result.g2$theta[-(1:500)]^(sum(data$delta)+a)*result.g2$tau[-(1:500)]^(sum(data$delta[data$group==1])+b)*exp(-result.g2$theta[-(1:500)]*(sum(data$time[data$group==2])+c)-result.g2$tau[-(1:500)]*result.g2$theta[-(1:500)]*(sum(data$delta[data$group==1])+d))
summary(post)
result.d2=matrix(ncol=4,nrow=2)
rownames(result.d2)=c("theta","tau")
colnames(result.d2)=c("mean","sd","95%lower","95%upper") 
result.d2[1,]=c(mean(result.g2$theta[-(1:500)]),sd(result.g2$theta[-(1:500)]),quantile(result.g2$theta[-(1:500)],c(0.25,0.75))) 
result.d2[2,]=c(mean(result.g2$tau[-(1:500)]),sd(result.g2$tau[-(1:500)]),quantile(result.g2$tau[-(1:500)],c( 0.25,0.75)))
##e
par(mfrow=c(1,1))
x=seq(0,6,by=0.01) 
plot(x,tau.prior(x,a,b,c,d)/integrate(tau.prior,0,Inf,a=a,b=b,c=c,d=d)$value,type="l",ylab="density",xlab=expression(tau),xlim=c(0,3),main="density of tau")
lines(density(result.g2$tau),lwd=2,lty=2,col=2)
legend(2,3,c("Prior","Posterior"),lwd=1:2,lty=1:2,col=1:2)
##f
t.test(result.g2$tau,mu=1,alternative="less",conf.level=0.95)
