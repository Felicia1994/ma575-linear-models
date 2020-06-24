
####### read data
day_data = read.csv('day.csv')
pairs(day_data)
par(mfrow=c(1,1))

####### select columns that we may use later, and transform windspeed
index = day_data$instant
season = day_data$season
year = day_data$yr
month = day_data$mnth
absmonth = abs(month-7)
holiday = day_data$holiday
weekday = day_data$weekday
workingday = day_data$workingday
weather = day_data$weathersit
temp = day_data$temp
adjtemp = day_data$atemp
absadjtemp = abs(adjtemp-0.665417)
# t_adjtemp = log(adjtemp)
humidity = day_data$hum
windspeed = day_data$windspeed
logwindspeed = log(windspeed)
sqrtwindspeed = sqrt(windspeed)
count = day_data$cnt

####### form useful columns into a new dataframe
data = data.frame(index,season,year,month,absmonth,holiday,weekday,workingday,weather,temp,adjtemp,absadjtemp,humidity,windspeed,logwindspeed,sqrtwindspeed,count)

#######
data_2011 = data[(1:365),]
data_2012 = data[(366:731),]
set.seed(1)
train = sample(1:nrow(data_2011), nrow(data_2011)/2)
train_2011 = data_2011[c(train),]
val = (-train)
val_2011 = data_2011[val,]

#######
pairs(data)

####### hist for each variable, giving the reason to transform a variable or not
hist(season)
hist(year)
hist(month)
hist(holiday)
hist(weekday)
hist(workingday)
hist(weather)
hist(temp)
hist(adjtemp)
hist(absadjtemp)
# hist(t_adjtemp)
hist(humidity)
par(mfrow=c(1,2))
hist(windspeed,breaks = seq(from = 0, to = 0.6, by = 0.01))
hist(logwindspeed,breaks = seq(from = -4, to = 0, by = 0.05))
hist(sqrtwindspeed,,breaks = seq(from = 0, to = 1, by = 0.01))
hist(count)

####### CoxBox method, giving why we use sqrt for windspeed
library(MASS)
Box = boxcox(windspeed~1,lambda=seq(-1,1,0.05))
Cox = data.frame(Box$x,Box$y)
Cox2 = Cox[with(Cox,order(-Cox$Box.y)),]
lambda = Cox2[1,'Box.x']
print(lambda)
T_box = (windspeed^lambda-1)/lambda
hist(T_box)

#Box = boxcox(adjtemp~count,lambda=seq(-1,1,0.05))
#Cox = data.frame(Box$x,Box$y)
#Cox2 = Cox[with(Cox,order(-Cox$Box.y)),]
#lambda = Cox2[1,'Box.x']
#print(lambda)
#T_box = (windspeed^lambda-1)/lambda
#hist(T_box)


####### simle linear model with sqrt transfromation for windspeed
# m3 = lm(count~season+year+month+I(month^2)+holiday+weekday+workingday+weather+adjtemp+I(adjtemp^2)+humidity+sqrtwindspeed,data = data)
m3 = lm(count~season+year+month+holiday+weekday+workingday+weather+adjtemp+humidity+sqrtwindspeed,data = data_2011)
m3 = lm(count~season+year+month+holiday+weekday+workingday+weather+absadjtemp+humidity+sqrtwindspeed,data = data_2011)
m3 = lm(count~season+absmonth+holiday+weekday+workingday+weather+absadjtemp+humidity+sqrtwindspeed,data = data_2011)
summary(m3)

#inverseResponsePlot(m3,lambda = c(0,1,0.1))
#invResPlot(m3)

library(car)
par(mfrow=c(1,1))
avPlot(m3,season)
#avPlot(m3,year)
#avPlot(m3,month)
avPlot(m3,absmonth)
avPlot(m3,holiday)
avPlot(m3,weekday)
avPlot(m3,workingday)
avPlot(m3,weather)
#avPlot(m3,adjtemp)
avPlot(m3,absadjtemp)
avPlot(m3,humidity)
avPlot(m3,sqrtwindspeed)

par(mfrow=c(2,2))
plot(m3)

####### test time correlation
StanRes3 <- rstandard(m3)
par(mfrow=c(1,1))
acf(StanRes3,main="Series Standardized Residuals")

####### gls
# install.packages('leaps')
library(leaps)
library(nlme)

# g3 <- gls(count~season+year+month+I(month^2)+holiday+weekday+workingday+weather+adjtemp+I(adjtemp^2)+humidity+sqrtwindspeed,correlation = corARMA(form = ~index),data=data,method="ML")
g3 <- gls(count~season+year+month+holiday+weekday+workingday+weather+adjtemp+humidity+sqrtwindspeed,correlation = corAR1(form = ~index),data=data_2011,method="ML")
g3 <- gls(count~season+year+month+holiday+weekday+workingday+weather+absadjtemp+humidity+sqrtwindspeed,correlation = corAR1(form = ~index),data=data_2011,method="ML")
g3 <- gls(count~season+absmonth+holiday+weekday+workingday+weather+absadjtemp+humidity+sqrtwindspeed,correlation = corAR1(form = ~index),data=data_2011,method="ML")
summary(g3)
# phi = g3$modelStruct$corStruct
phi = 0.5447225
phi = 0.365485
phi = 0.3893127
phi = 0.3769957

m3 = lm(count~season+absmonth+holiday+weekday+workingday+weather+absadjtemp+humidity+sqrtwindspeed,data = data_2011)
x <- model.matrix(m3)
Sigma <- diag(365)
Sigma <- phi^abs(row(Sigma)-col(Sigma))
sm <- chol(Sigma)
smi <- solve(t(sm))
xstar <- smi %*% x
ystar <- smi %*% count[c(1:365)]
m1tls <- lm(ystar ~ xstar-1) 
summary(m1tls)
par(mfrow=c(2,2))
plot(m1tls)

####### no time correlation after gls
StanRes3 <- rstandard(m1tls)
par(mfrow=c(1,1))
acf(StanRes3,main="Series Standardized Residuals")

####### train
m3 = lm(count~season+absmonth+holiday+weekday+workingday+weather+absadjtemp+humidity+sqrtwindspeed,data = train_2011)
x <- model.matrix(m3)
Sigma <- diag(length(count[train]))
Sigma <- phi^abs(row(Sigma)-col(Sigma))
sm <- chol(Sigma)
smi <- solve(t(sm))
xstar <- smi %*% x
ystar <- smi %*% count[train]
m1tls <- lm(ystar ~ xstar-1) 
summary(m1tls)
par(mfrow=c(2,2))
plot(m1tls)

#######
X <- xstar
Y <- ystar
b <- regsubsets(as.matrix(X), Y) 
rs <- summary(b)
print(rs)

###### val
m3 = lm(count~season+absmonth+holiday+weekday+workingday+weather+absadjtemp+humidity+sqrtwindspeed,data = val_2011)
x <- model.matrix(m3)
Sigma <- diag(365-length(count[train]))
Sigma <- phi^abs(row(Sigma)-col(Sigma))
sm <- chol(Sigma)
smi <- solve(t(sm))
xstar <- smi %*% x
ystar <- smi %*% count[c(1:365)][val]
m1tls <- lm(ystar ~ xstar-1) 
summary(m1tls)
par(mfrow=c(2,2))
plot(m1tls)

X <- xstar
Y <- ystar
b <- regsubsets(as.matrix(X), Y) 
rs <- summary(b)
print(rs)

om1<- lm(count~ X[,8],data=val_2011)
om2<- lm(count~ X[,8]+X[,7],data=val_2011)
om3<- lm(count~ X[,8]+X[,7]+X[,2],data=val_2011)
om4<- lm(count~ X[,8]+X[,7]+X[,2]+X[,10],data=val_2011)
om5<- lm(count~ X[,8]+X[,7]+X[,2]+X[,10]+X[,1],data=val_2011)
om6<- lm(count~ X[,8]+X[,7]+X[,2]+X[,10]+X[,1]+X[,4],data=val_2011)
om7<- lm(count~ X[,8]+X[,7]+X[,2]+X[,10]+X[,1]+X[,4]+X[,9],data=val_2011)
om8<- lm(count~ X[,8]+X[,7]+X[,2]+X[,10]+X[,1]+X[,4]+X[,9]+X[,5],data=val_2011)

n <- nrow(data)
lm.stat <- function(l) {
  npar = length(coef(l)) + 1 
  AIC = extractAIC(l, k = 2)[2] 
  BIC = extractAIC(l, k = log(n))[2] 
  return( c(rs$adjr2[npar - 2], #adjr2 
            AIC, 
            AIC + 2 * npar * (npar + 1) / (n - npar + 1), # AICc 
            BIC))
}
matrix(unlist(lapply(list(om1, om2, om3, om4, om5, om6, om7,om8),lm.stat)), byrow = TRUE, ncol = 4, dimnames = list(1:8, c("adjr2", "AIC", "AICc", "BIC")))

par(mfrow=c(2,2))
plot(om4)
###### prediction
m3 = lm(count~season+absmonth+holiday+weekday+workingday+weather+absadjtemp+humidity+sqrtwindspeed,data = data_2012)
x <- model.matrix(m3)
Sigma <- diag(366)
Sigma <- phi^abs(row(Sigma)-col(Sigma))
sm <- chol(Sigma)
smi <- solve(t(sm))
xstar <- smi %*% x

X <- xstar

prediction_raw = predict(om4,data.frame(X[,8]+X[,7]+X[,2]+X[,9]+X[,10]))
count_pred_raw = solve(smi) %*% prediction_raw
sqrt(mean((count[(366:731)]-count_pred_raw)^2))

######
avg_1 = mean(count[data[4]==1&data[3]==0])
avg_2 = mean(count[data[4]==2&data[3]==0])
avg_3 = mean(count[data[4]==3&data[3]==0])
avg_4 = mean(count[data[4]==4&data[3]==0])
avg_5 = mean(count[data[4]==5&data[3]==0])
avg_6 = mean(count[data[4]==6&data[3]==0])
avg_7 = mean(count[data[4]==7&data[3]==0])
avg_8 = mean(count[data[4]==8&data[3]==0])
avg_9 = mean(count[data[4]==9&data[3]==0])
avg_10 = mean(count[data[4]==10&data[3]==0])
avg_11 = mean(count[data[4]==11&data[3]==0])
avg_12 = mean(count[data[4]==12&data[3]==0])

avg_1 
avg_2 
avg_3 
avg_4 
avg_5 
avg_6 
avg_7 
avg_8 
avg_9 
avg_10
avg_11
avg_12


count_pred = count_pred_raw*0.93+733
sqrt(mean((count[(366:731)]-count_pred)^2))

avg_12/avg_1

par(mfrow=c(1,1))
plot(count_pred_raw,count[(366:731)],pch=20)
lm(count[(366:731)]~count_pred_raw)

plot(index[(366:731)],count[(366:731)],col='blue',pch=20)
points(index[(366:731)],count_pred_raw,col='red',pch=20)
legend(500,2000, legend=c('real count','predicted count'), col=c('blue','red'), lty=0, cex=1, pch=1)



