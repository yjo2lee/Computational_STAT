setwd("~/Desktop/Grad2/Computational_Stat/HW/datasets/")
library(glmulti)
library(metafor)
library(caret)
library(subselect)
library(GenSA)
library(gensa)
library(GA)
library(Metrics)

data = read.table("baseball.dat", header = TRUE)
# predictors = subset(data,select=-salary)
# data["logsalary"] = log(data$salary)
# test_data = data[2:29,]

# simulated annealing
## GenSA 이용
SA.ftn <- function(par){
  par.var = data[,-1][,par>0.5]
  lm = lm(log(data$salary)~., data=par.var)
  AIC = extractAIC(lm)[2]
  return(AIC)
}

#sa1.model = GenSA(fn=SA.ftn, lower=rep(0, 12), upper=rep(1, 12), control=list(max.call=2000, smooth=F, verbose=T))
sa2.model = GenSA(fn=SA.ftn, lower=rep(0, 27), upper=rep(1, 27), control=list(max.call=2000, smooth=F, verbose=T))
test_lm2 = lm(log(data$salary) ~ 1 + obp + homeruns + walks + sbs + errors + freeagent + arbitration + runsperso + soserrors + sbsobp+sbshits , data = data)
system.time(GenSA(fn=SA.ftn, lower=rep(0, 27), upper=rep(1, 27), control=list(max.call=2000, smooth=F, verbose=T))) 

# rmse
train.size = as.integer(length(data$salary)*0.7)
train.data = data[1:train.size,]
test.data = data[train.size:length(data$salary),]
model.1 = lm(log(salary) ~ obp + homeruns + walks + sbs + errors + freeagent + arbitration + runsperso + soserrors + sbsobp + sbshits, data = train.data)
model.2 = lm(log(salary) ~ average + obp + doubles + homeruns + 
               walks + errors + freeagent + arbitration + runsperso + hrspererror + 
               soserrors + sbsobp, data = train.data)
predictions.1 = predict(model.1, test.data)
predictions.2 = predict(model.2, test.data)
rmse.1 = rmse(log(test.data$salary), predictions.1)
rmse.2 = rmse(log(test.data$salary), predictions.2)
## anneal 이용
# logmat = lmHmat(logsalary ~ ., data = test_data)
# fullmat = lmHmat(salary ~., data = data[1:28, ])
# sa1_results = anneal(logmat$mat, kmin = 10, kmax = 13, nsol =3, niter=10, H=logmat$H, r=logmat$r, criterion="RM")
# sa2_results = anneal(fullmat$mat, kmin = 10, kmax = 13, nsol =3, niter=10,criterion="RM")
print(sa1_results$bestsets)


### genetic algorithm
x_col = colnames(data)[2:length(colnames(data))]
g1 <- glmulti("salary", xr =  x_col, data=data,level=1, crit = 'aicc')
g2 <- glmulti("logsalary", xr=x_col, data=data, level=1, method = 'g', crit = 'aicc', marginality=TRUE)
g3 <- glmulti("logsalary", xr=x_col, data=data, level=1, method = 'g', crit = 'aic', marginality=FALSE, maxsize = 12, popsize = 20, mutrate = 0.1)

GA.ftn <- function(run){
  run.vars <- data[,-1][, run == 1]
  lm = lm(log(data$salary)~., data=run.vars)
  AIC = extractAIC(lm)[2]
  return(-AIC)
}
m = 27
GA <- ga(type="binary", fitness=GA.ftn, lower=rep(0, m), upper=rep(1, m), nBits=m)
system.time(ga(type="binary", fitness=GA.ftn, lower=rep(0, m), upper=rep(1, m), nBits=m)) 
summary(GA)


mykmeans<-function(x,k,nstart=1){
  totwithin.old = Inf
  for(n in 1:nstart){
    flag = TRUE 
    in.mean = x[sample(nrow(x),k),] 
    member = rep(NA,nrow(x)) 
    cluster.mean = in.mean
    #find center 
    while(flag){
      for(i in 1:nrow(x)){
        #거리계산 
        distance1<-as.matrix(dist(rbind(x[i,],in.mean)))[1,-1] 
        #membership update 
        member[i]<-which.min(distance1)
      }
      for(i in 1:k){cluster.mean[i,]<-colMeans(x[member==i,]) } 
      ifelse(sum(in.mean-cluster.mean)==0,flag<-FALSE,in.mean <-cluster.mean)
    }
    #minimize TWSS 
    withinss<-rep(0,k)
    for(i in 1:k){
      withinss[i]<-sum(as.matrix(dist(rbind(x[member==i,],cluster.mean[i,]))^2)[nrow(x[member==i,])+1,- nrow(x[member==i,])-1])
    } 
    totwithinss<-sum(withinss) 
    if(totwithinss<totwithin.old){
      cluster.list<-list(cluster.membership=member,cluster.mean=cluster.mean,Withinss=withinss,TotwithinSS=totwithinss)
      totwithin.old<-totwithinss }
  }
  return(cluster.list) 
}


test1 = rep(0, 30)
test1[1:10]  = rnorm(10, 0, 1)
test1[11:20] = rnorm(10, 3, 2)
test1[21:30] = rnorm(10, 5, 3)
test2 = rep(0, 30)
test2[0:20]  = rnorm(10, 0, 1)
test2[11:20] = rnorm(10, 3, 2)
test2[21:30] = rnorm(10, 5, 3)
test_data = data.frame(test1 = test1, test2 = test2)
#test_data = data.frame(test = test)
mykmeans(test_data, 3, nstart = 10)
kmeans(test_data, 3, nstart = 10)
test1
system.time(mykmeans(test_data,3,nstart=10)) 
system.time(for(i in 1:100){kmeans(test_data,2,nstart=10)})/100
