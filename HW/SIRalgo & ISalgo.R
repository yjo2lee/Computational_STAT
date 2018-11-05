# slash distribution

X <- dnorm(0, 1)
U <- dunif(0, 1)

slash_density<-function(y){
  if (y == 0){
    d = 1/(2*sqrt(2*pi))
  }
  else {
    d = f.y <- (1-exp((-y^2)/2))/(y^2*sqrt(2*pi)
  }
  return(d)
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