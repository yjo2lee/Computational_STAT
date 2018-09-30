setwd("~/Desktop/Grad2/Computational_Stat/HW/datasets/")
library(MASS)
library(leaps)
library(Metrics)
data = read.table("baseball.dat", header = TRUE)

#1) Stepwise
fullmodel = lm(salary~., data = data)
stepwise.both = step(fullmodel, direction = "both")
stepwise.backward = step(fullmodel, direction = "backward")
stepwise.forward = step(fullmodel, direction = "forward")

AIC(stepwise.both)
AIC(stepwise.backward)
AIC(stepwise.forward)

plot(stepwise.both)
# stepwise.2 = stepAIC(fullmodel, direction = "both")
# stepwise.3 = extractAIC(lm(salary ~ runs + hits + rbis + sos + sbs + freeagent + 
#                          arbitration + runsperso + hitsperso + hrsperso + rbisperso + 
#                          walksperso + soserrors + sbsobp, data = data))

# log(y): transformation
fullmodel.trans = lm(log(salary)~., data = data)
stepwise.trans.both = step(fullmodel.trans, direction = "both")
AIC(stepwise.trans.both)
plot(stepwise.trans.both)
# stepwise.trans.2 = stepAIC(fullmodel.trans, direction = "both")
# stepwise.trans.3 = extractAIC(lm(log(salary) ~ average + runs + triples + rbis + 
#                 sos + freeagent + arbitration + runsperso + hitsperso + soserrors + 
#                 sbsobp + sbsruns, data = data))

#2) All possible regression
library(leaps)
allposs = regsubsets(salary~., data = data, nbest = 1, nvmax = NULL, method = "exhaustive")
summary(allposs)$cp

AIC_all_poss = rep(0,27)
salary = data$salary
lm_model = list(NA)
for(i in 1:27){
  lm_model[[i]] = lm(salary~., data = data[summary(allposs)$which[i,]])
  AIC_all_poss[i] = AIC(lm_model[[i]])
}
lm_model[[which.min(AIC_all_poss)]]

# all.possible.regression_log(salary)
allposs_trans = regsubsets(log(salary)~., data = data, nbest = 1, nvmax = NULL, method = "exhaustive")
AIC_all_poss_trans = rep(0,27)
lm_model_trans = list(NA)
for(i in 1:27){
  lm_model_trans[[i]] = lm(log(salary)~., data = data[summary(allposs_trans)$which[i,]])
  AIC_all_poss_trans[i] = AIC(lm_model_trans[[i]])
}
allposs.best = lm_model_trans[[which.min(AIC_all_poss_trans)]]

# stepwise vs all.possible.reg
system.time(for(i in 1:10) step(lm(log(salary)~., data = data), direction = "both"))
system.time(for(i in 1:10) regsubsets(log(salary)~., data = data, nbest = 1, nvmax = NULL, method = "exhaustive"))

# rmse
train.size = as.integer(length(data$salary)*0.7)
train.data = data[1:train.size,]
test.data = data[train.size:length(data$salary),]
model.1 = lm(log(salary) ~ average + runs + triples + rbis + sos + freeagent + arbitration + runsperso + hitsperso + soserrors + sbsobp + sbsruns, data = train.data)
model.2 = lm(log(salary) ~ obp + runs + triples + rbis + sos + freeagent + arbitration + runsperso + hitsperso + soserrors + sbsobp + sbsruns, data = train.data)
predictions.1 = predict(model.1, test.data)
predictions.2 = predict(model.2, test.data)
rmse.1 = rmse(log(test.data$salary), predictions.1)
rmse.2 = rmse(log(test.data$salary), predictions.2)