# multinomial with complex cell structure
multiEM_moth<-function(nc,ni,nt,threshold=10^(-10)){ 
  p = 1/3 ; q = 1/3 ; r = 1-p-q; n = nc+ni+nt
  error1 = 1 ; error2 = 1; maxiter = 1000; iter = 0
  while(error1>=threshold && error2>=threshold&& iter < maxiter){ 
    #E-Step
    ncc = nc*p^2/(p^2+2*p*r+2*p*q) 
    nci = nc*2*p*q/(p^2+2*p*r+2*p*q) 
    nct = nc*2*p*r/(p^2+2*p*r+2*p*q) 
    nii = ni*q^2/(q^2+2*q*r) 
    nit = ni*2*q*r/(q^2+2*q*r)
    
    nc.p = ncc+0.5*nci+0.5*nct 
    ni.p = nii+0.5*nci+0.5*nit 
    nr.p = nt+0.5*nit+0.5*nci
    #M-Step
    p.0 = p ; q.0<-q
    p = nc.p/n ; q = ni.p/n ; r = 1-p-q
    #Stopping Rule
    error1 = abs(p-p.0) ; error2 = abs(q-q.0); iter = iter + 1
  }
  return(list(iteration=iter,prob=c(p=p,q=q,r=r),freq=c(CC=ncc,CI=nci,CT=nct,II=nii,IT=nit, TT=nt)))
}
multiEM_moth(85,196,341)
system.time(for(i in 1:10^3) multiEM_moth(85,196,341))

# Newton cotes
# Riemann rule
g1<- function(x){
  return(x^2)
}
g2 <- function(x){
  cos(x+1/x)/sqrt(x)
}
g3 <- function(x){
  ((3^(-x)/(x^2+6))+x^3)*sin(x)
}

myRiemann<- function(f, a, b, threshold = 10^(-5)){
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
round(myRiemann(f1, 0,2),4)

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

plot(seq(0,5,0.1),g1(seq(0,5,0.1)),type='l',xlab="x",ylab="g1(x)",main="g1")
round(myRiemann(g1, 0, 2),4)
round(myTPz(g1, 0, 2),4)
round(mySimpson(g1, 0, 2),4)
integrate(g1,0,2)

plot(seq(1,10,0.1),g2(seq(1,10,0.1)),type='l',xlab="x",ylab="g2(x)",main="g2")
round(myRiemann(g2, 2,6),4)
round(myTPz(g2, 2,6),4)
round(mySimpson(g2, 2,6),4)
integrate(g2,2,6)

plot(seq(1,10,0.1),g3(seq(1,10,0.1)),type='l',xlab="x",ylab="g3(x)",main="g3")
round(myRiemann(g3, 4, 7),4)
round(myTPz(g3, 4, 7),4)
round(mySimpson(g3, 4, 7),4)
integrate(g3,4,7)

system.time(myRiemann(g3, 4, 7))
system.time(for(i in 1:10^3) myTPz(g3, 4, 7))
system.time(for(i in 1:10^3) mySimpson(g3, 4, 7))
system.time(for(i in 1:10^3) integrate(g3,4,7))

