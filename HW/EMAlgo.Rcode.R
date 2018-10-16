# univariate 2-class GMM
myGMM <- function(x, threshold = 10^(-10)){
  ## initialize
  pi0 = 1/2
  mu = sample(x,2); mu1 = min(mu); mu2 = max(mu)
  sigma1 = sd(x); sigma2 = sd(x)
  iter = 0; maxiter = 30000; error = 1
  while(error >= threshold && iter < maxiter){
    ## E-step: Y update
    y = dnorm(x,mu1,sigma1)*pi0/(dnorm(x,mu1,sigma1)*pi0+dnorm(x,mu2,sigma2)*(1-pi0))
    ## M-step: theta update
    pi1 = mean(y)
    mu1 = sum(y*x)/sum(y); mu2 = sum((1-y)*x)/sum(1-y)
    sigma1 = sqrt(sum(y*(x-mu1)^2)/sum(y)); sigma2 = sqrt(sum((1-y)*(x-mu2)^2)/sum(1-y))
    error = abs(pi0 - pi1)
    iter = iter + 1
    pi0 = pi1
  }
  y1 = length(y[y>0.5]); y2 = length(y[y<= 0.5])
  return (data.frame(y1,y2, pi1, mu1, mu2, sigma1, sigma2, error, iter))
}
## part1_modelA
x1 = rnorm(100, 0, 1); x2 = rnorm(400, 5, 3); x = c(x1, x2)
x1 = rnorm(1000, 0, 1); x2 = rnorm(4000, 5, 3); x = c(x1, x2)
myGMM(x)
system.time(myGMM(x))
plot(density(x1), xlim = c(-5, 12), ylim = c(0,0.4),xlab = 'x', ylab = 'density', main = 'density of N(0,1), N(5,3), GMM_modelB')
lines(density(x2))
lines(density(x), col = 'red')

## part1_modelC
x1 = rnorm(1000, 0, 2); x2 = rnorm(4000, 2, 4); x = c(x1, x2)
myGMM(x)
system.time(myGMM(x))
plot(density(x1), xlim = c(-5, 12), ylim = c(0,0.4),xlab = 'x', ylab = 'density', main = 'density of N(2,2), N(4,4), GMM_modelC')
lines(density(x2))
lines(density(x), col = 'red')

# Bivariate normal with missing values
w1 = c(8,11,16,18,6,4,20,25,9,13) 
w2 = c(10,14,16,15,20,4,18,22,NA,NA)
W = data.frame(w1, w2)

mybivariate <- function(w1, w2, threshold = 10^(-10)){
  iter = 0; maxiter = 1000; error = 1; mu2.0 = mean(w2[1:8]); cov = cov(w1[1:8], w2[1:8])
  mu1 = mean(w1); sigma1 = (sum(w1^2)-sum(w1)^2/length(w1))/length(w1)
  while(error >= threshold && iter < maxiter){
    #E-step
    w2[9] =  mu2.0 + (cov/sigma1)*(w1[9]-mu1)
    w2[10] = mu2.0 + (cov/sigma1)*(w1[10]-mu1)
    
    # M-step
    mu2.1 = mean(w2)
    sigma2 = (sum(w2^2)-sum(w2)^2/length(w2))/length(w2)
    cov = (sum(w1*w2)-sum(w1)*sum(w2)/length(w1))/length(w1)
    error = abs(mu2.1-mu2.0)
    mu2.0 = mu2.1
    iter = iter + 1
  }
  return(data.frame(w2[9], w2[10], mu1, mu2.1, sigma1, sigma2, cov, error, iter))
}
View(round(mybivariate(w1,w2),4))
system.time(mybivariate(w1,w2))

# regression modelling with missing values
setwd("~/Desktop/Grad2/Computational_Stat/HW/")
data3 = read.table("auto_mpg_data.txt", sep = "")
data3 = data3[,-c(4,9)]
colnames(data3) = c("mpg","cylinders","displacement","weight","acceleration","m odel_year","origin")
data3$origin = as.factor(data3$origin)
data3$cylinders = as.factor(data3$cylinders)

missing.percent = seq(0.01,0.1, 0.01)
n = length(data3$mpg)
del_index = sample(n,as.integer(n*missing.percent[10]))
raw_data3 = data3
trueAIC = AIC(lm(mpg~., raw_data3))
true.beta = lm(mpg~., raw_data3)$coef
data3[del_index,1] = NA

missingreg <- function(data, del_index, threshold=10^(-10)){
  data[del_index,1] = mean(data$mpg, na.rm = TRUE)
  beta.0 = lm(mpg~., data)$coef
  iter = 0; maxiter = 1000; error = 1; 
  while(error >= threshold && iter < maxiter){
    data[del_index,1] = lm(mpg~., data)$fitted.values[del_index]
    AIC = AIC(lm(mpg~., data))
    beta.1 = lm(mpg~., data)$coef
    #print(data$mpg[del_index])
    fity = predict(lm(mpg~., data))
    #fity = lm(mpg~., data)$fitted.values
    #print(fity)
    error = sum(abs(beta.1-beta.0)^2)
    beta.0 = beta.1
    iter = iter + 1
  }
  return(list(fity, beta.1, error, iter, AIC))
}

fit_list = missingreg(data3, del_index)
fity = fit_list[[1]]
sel.y = fity[del_index]
truey = lm(mpg~.,raw_data3)$fitted.values[del_index]
View(round(data.frame(sel.y, truey),4))

embeta = fit_list[[2]]
View(round(data.frame(embeta, true.beta),2))

system.time(missingreg(data3, del_index)) 
plot(sel.y,truey, xlab="yhat",ylab="truey") 
abline(0,1,col=2)

# multinomial with complex cell structure
multiEM<-function(no,na,nb,nab,threshold=10^(-10)){ 
  p = 1/3 ; q = 1/3 ; r = 1-p-q; n = no+na+nb+nab
  error1 = 1 ; error2 = 1; maxiter = 1000; iter = 0
  while(error1>=threshold && error2>=threshold&& iter < maxiter){ 
    #E-Step
    naa = na*p^2/(p^2+2*p*r) 
    nao = na*2*p*r/(p^2+2*p*r) 
    nbb = nb*q^2/(q^2+2*q*r) 
    nbo = nb*2*q*r/(q^2+2*q*r) 
    na.p = naa+0.5*nao+0.5*nab 
    nb.p = nbb+0.5*nbo+0.5*nab 
    no.p = no+0.5*nao+0.5*nbo
    #M-Step
    p.0 = p ; q.0<-q
    p = na.p/n ; q = nb.p/n ; r = 1-p-q
    #Stopping Rule
    error1 = abs(p-p.0) ; error2 = abs(q-q.0); iter = iter + 1
  }
  return(list(iteration=iter,prob=c(p=p,q=q,r=r),freq=c(AA=naa,AO=nao,BB=nbb,BO=nbo,AB=nab, OO=no)))
}
multiEM(176,182,60,17)
system.time(for(i in 1:10^3) multiEM(176,182,60,17))





