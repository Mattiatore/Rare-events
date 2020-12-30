# loading data
library(ismev);
library(evd);
load(file="Project2020.Rdata")

# allow to print entire dataset
options(max.print=2000)

# clear seasonality, maxima are during summer
par(mfrow=c(1,1))
plot(ts(temp_beznau$max))
savePlot("maxima.jpg")

#remove 2013 since only two months?
temp_beznau=temp_beznau[-2:-1,]

# extract annual maxima and their return level
annual<-integer((nrow(temp_beznau)-11)/12+1)
annual[1]<-max(temp_beznau$max[1:11])
annual[2:length(annual)]<-sapply(split(temp_beznau$max[12:nrow(temp_beznau)], rep(1:(length(temp_beznau$max[12:nrow(temp_beznau)])/12), each=12)), max)

#fit for annual maxima
annual_fit = gev.fit(annual)

# estimated errors for return levels
return_statio_error <- function(r){
  y = -log(1-(1/r))
  del=matrix(ncol=1,nrow=3)
  del[1,1]=1
  del[2,1]=-((annual_fit$mle[3])^(-1))*(1-(y^(-annual_fit$mle[3])))
  del[3,1]=((annual_fit$mle[2])*((annual_fit$mle[3])^(-2))*(1-((y)^(-annual_fit$mle[3]))))-((annual_fit$mle[2])*((annual_fit$mle[3])^(-1))*((y)^(-(annual_fit$mle[3])))*log(y))
  del.transpose=t(del)
  return (sqrt(del.transpose%*%annual_fit$cov%*%del))
}

#r-year return level stationary case
return_statio <- function(r){
  return (annual_fit$mle[1]+annual_fit$mle[2]/annual_fit$mle[3]*(-1+(-log(1-1/r))**(-annual_fit$mle[3])))
}
print("Return levels for annual maxima based GEV and estimated standar error.")
print(paste("100 year return level:",toString(return_statio(100)),"estiamted standard error",toString(return_statio_error(100)),sep=" "))
print(paste("1000 year return level:",toString(return_statio(1000)),"estiamted standard error",toString(return_statio_error(1000)),sep=" "))
print(paste("10000 year return level:",toString(return_statio(10000)),"estiamted standard error",toString(return_statio_error(10000)),sep=" "))

# fitting GEV for stationary case
monthly_fit = gev.fit(temp_beznau$max)
n <- nrow(temp_beznau)
AIC_mon <- 2*monthly_fit$nllh+2*(2*K+3) +2*(2*K+3)*(2*K+4)/(n-2*K-4)

# model with trig functions
# K from 1 to 6
# assign t0 as in the description 
t0 = as.Date("2000/01/01")
diff <- seq(nrow(temp_beznau))
for (i in seq(nrow(temp_beznau))) {
  str <- paste(toString(temp_beznau$year[i]),toString(temp_beznau$mon[i]),toString(temp_beznau$day[i]),sep="/")
  diff[i] <- as.Date(str)-t0
}
AIC_monthly <- seq(6)
for (K in seq(6)){
  t <- matrix(ncol=K*2+1,nrow=nrow(temp_beznau))
  t[,1] <- diff/(100*365.25)
  # sine 
  for (i in seq(K)) {
    t[,i+1] <- sin(2*i*pi*diff/365.25)
  }
  # cos
  for (i in seq(K)) {
    t[,i+K+1] <- cos(2*i*pi*diff/365.25)
  }
  monthly_trig <- gev.fit(temp_beznau$max,ydat=t,mul=c(1:(2*K+1)))
  
  # AIC = -2 * LL + 2 * p + 2*p*(p+1)/(n-p-1)
  AIC_monthly[K] <- 2*monthly_trig$nllh+2*(2*K+3) +2*(2*K+3)*(2*K+4)/(n-2*K-4)
}
t <- matrix(ncol=1,nrow=nrow(temp_beznau))
t[,1] <- diff/(100*365.25)
monthly_trig <- gev.fit(temp_beznau$max,ydat=t,mul=c(1))
K <- 0
AIC_monthly1 <- 2*monthly_trig$nllh+2*(2*K+3) +2*(2*K+3)*(2*K+4)/(n-2*K-4)

# plot diagnostic for constant
gev.diag(monthly_trig)

# plot for best K=1
K <- 1
t <- matrix(ncol=K*2+1,nrow=nrow(temp_beznau))
t[,1] <- diff/(100*365.25)
# sine 
for (i in seq(K)) {
  t[,i+1] <- sin(2*i*pi*diff/365.25)
}
# cos
for (i in seq(K)) {
  t[,i+K+1] <- cos(2*i*pi*diff/365.25)
}
monthly_trig <- gev.fit(temp_beznau$max,ydat=t,mul=c(1:(2*K+1)))

# plot diagnostic for K=6 from (1)
gev.diag(monthly_trig)

if((monthly_trig$mle[1]-1.96*monthly_trig$se[1])<0){
  print("There is no clear trend because u_1 is not different from 0 at the 95% CI")
} else{
  print("There is a trend because u_1 is different from 0 at the 95% CI")
}

# there is seasonal variation because some coefficients of the sin and cos are different from 0 at the 95% CI
# test model without linear term u_1 in equation (1)
K <- 1
t <- matrix(ncol=K*2,nrow=nrow(temp_beznau))
# sine 
for (i in seq(K)) {
  t[,i] <- sin(2*i*pi*diff/365.25)
}
# cos
for (i in seq(K)) {
  t[,i+K] <- cos(2*i*pi*diff/365.25)
}
monthly_trig <- gev.fit(temp_beznau$max,ydat=t,mul=c(1:(2*K)))
AIC_month <- 2*monthly_trig$nllh+2*(2*K) +2*(2*K)*(2*K+1)/(n-2*K-1)

# plot diagnostic for best model without U_1
gev.diag(monthly_trig)
savePlot("diag_best.jpg")

# return levels for non stationary case
return_non_statio <- function(r,t){
  L=length(monthly_trig$mle)
  T <- as.numeric((as.Date(t)-t0), unit="days")
  somma <- monthly_trig$mle[2]*sin(2*pi*(T)/365.25) + monthly_trig$mle[3]*cos(2*pi*(T)/365.25)
  return (somma + monthly_trig$mle[1]+monthly_trig$mle[L-1]/monthly_trig$mle[L]*(-1+(-log(1-1/r))**(-monthly_trig$mle[L])))
}

# approximate standard error considering parameters stationary
return_error <- function(r){
  y = -log(1-(1/r))
  del=matrix(ncol=1,nrow=3)
  del[1,1]=1
  del[2,1]=-((monthly_trig$mle[5])^(-1))*(1-(y^(-monthly_trig$mle[5])))
  del[3,1]=((monthly_trig$mle[4])*((monthly_trig$mle[5])^(-2))*(1-((y)^(-monthly_trig$mle[5]))))-((monthly_trig$mle[4])*((monthly_trig$mle[5])^(-1))*((y)^(-(monthly_trig$mle[5])))*log(y))
  del.transpose=t(del)
  return (sqrt(del.transpose%*%monthly_trig$cov[-3:-4,-3:-4]%*%del))
}

# print errors
print("Return levels for annual maxima based on non-stationary GEV for day 15 of each month in year 2030.")
for (i in seq(12)){
  print(paste("Month",toString(i),sep=" "))
  print(paste("100 year return level:",toString(return_non_statio(100,paste("2030",toString(i),"15",sep="/"))),"estiamted standard error",toString(return_error(100)),sep=" "))
  print(paste("1000 year return level:",toString(return_non_statio(1000,paste("2030",toString(i),"15",sep="/"))),"estiamted standard error",toString(return_error(1000)),sep=" "))
  print(paste("10000 year return level:",toString(return_non_statio(10000,paste("2030",toString(i),"15",sep="/"))),"estiamted standard error",toString(return_error(10000)),sep=" "))
}

# profile log-likelihoods for shape parameters of best model K=1 without u_1
gev.profxi(monthly_fit,-0.6,-0.2,conf=0.95,nint=100)
savePlot("profile_shape.jpg")

# Modelling pairwise dependence in the extremes
# drop last two as previously
temp_muehleberg=temp_muehleberg[-2:-1,]
# comparing with station Muhleberg
par(mfrow=c(1,1))
plot(ts(temp_muehleberg$max))
savePlot("maxima_mue.jpg")

# apply previous best model
K <- 1
t <- matrix(ncol=K*2,nrow=nrow(temp_muehleberg))
# change diff variable
diff <- seq(nrow(temp_muehleberg))
for (i in seq(nrow(temp_muehleberg))) {
  str <- paste(toString(temp_muehleberg$year[i]),toString(temp_muehleberg$mon[i]),toString(temp_muehleberg$day[i]),sep="/")
  diff[i] <- as.Date(str)-t0
}
# sine 
for (i in seq(K)) {
  t[,i] <- sin(2*i*pi*diff/365.25)
}
# cos
for (i in seq(K)) {
  t[,i+K] <- cos(2*i*pi*diff/365.25)
}
monthly_mue <- gev.fit(temp_muehleberg$max,ydat=t,mul=c(1:(2*K)))

if(monthly_mue$mle[5]+1.96*monthly_mue$se[5]>monthly_trig$mle[5]-1.96*monthly_trig$se[5]){
  print("The 95% CI for xi overlap in the two series, therefore we can assume the same xi")
} else{
  print("The 95% CI for xi do not overlap in the two series, therefore we cannot assume the same xi")
}

#put data togheter
temp <- matrix(ncol=2, nrow=nrow(temp_beznau))
temp[,1] <- temp_beznau$max
temp[,2] <- temp_muehleberg$max

# We now examine the asymptotic dependance by using chi-chibar plots
par(mfrow=c(1,2))
label_chi = sprintf("Chi plot")
label_chibar = sprintf("Chi Bar plot")
chiplot(temp, main1 = label_chi, which=1)
abline(a=0,b=0,col="red")
chiplot(temp, main2 = label_chibar, which=2)
abline(a=0,b=0,col="red")


#high correlation
cor(temp[,1],temp[,2])

# Fitting a bivariate extreme value model with logistic dependence structure
fit1 <- fbvevd(temp,model="log") 
fit1 
par(mfrow=c(3,2)) 
plot(fit1)
aic1 <- fit1$dev + 2*length(fit1$param)

#negative logistic
fit2 <- fbvevd(temp,model="neglog") 
fit2
par(mfrow=c(3,2)) 
plot(fit2)
aic2 <- fit2$dev + 2*length(fit2$param)

#bivariate Coles-Tawn 
fit3 <- fbvevd(temp,model="ct") 
fit3
par(mfrow=c(3,2)) 
plot(fit3)
aic3 <- fit3$dev + 2*length(fit3$param)

#bivariate negative logistic 
fit4 <- fbvevd(temp,model="negbilog") 
fit4
par(mfrow=c(3,2)) 
plot(fit4)
aic4 <- fit4$dev + 2*length(fit4$param)

# Following the AIC metric the "best" model on this dataset seems to be Coles-Tawn
print(fit3$dep.summary)
# alpha value seem to imply extremal dependence 

#sharing xi parameter
fit5 <- fbvevd(temp,model="ct",cshape=TRUE)
aic5 <- fit5$dev + 2*length(fit5$param)

# Based on that, we better have to use a common shape parameter (the two models cannot be considered significantly different).
fit5$dev-fit3$dev # -16.45
qchisq(0.95,1) # 3.84

# alternative, fit gev to each month separately
# return level considering only one month maxima series
return_level = function(x,period=100){
  if (is.list(x)) {
    loc = x$estimate[[1]]
    scale = x$estimate[[2]]
    shape = x$estimate[[3]]
  }
  if (is.vector(x)) {
    loc = x[1]
    scale = x[2]
    shape = x[3]
  }
  p = 1/period
  
  level = loc + scale*(((-log(1-p))^-shape-1)/shape)
  return(level)
}

#january
print("January:")
jan <- seq(10,311,12)
jan_fit <- fgev(temp_beznau$max[jan], 
                method = "Nelder-Mead")
par(mfrow=c(2,2))
plot(jan_fit)
print("Return level 100:")
return_level(jan_fit,100)
print("Return level 1000:")
return_level(jan_fit,1000)
print("Return level 10000:")
return_level(jan_fit,10000)

#february
print("February:")
jan <- seq(11,311,12)
jan_fit <- fgev(temp_beznau$max[jan], 
                method = "Nelder-Mead")
par(mfrow=c(2,2))
plot(jan_fit)
print("Return level 100:")
return_level(jan_fit,100)
print("Return level 1000:")
return_level(jan_fit,1000)
print("Return level 10000:")
return_level(jan_fit,10000)

#march
print("March:")
jan <- seq(12,311,12)
jan_fit <- fgev(temp_beznau$max[jan], 
                method = "Nelder-Mead")
par(mfrow=c(2,2))
plot(jan_fit)
print("Return level 100:")
return_level(jan_fit,100)
print("Return level 1000:")
return_level(jan_fit,1000)
print("Return level 10000:")
return_level(jan_fit,10000)

#other months
names <- c("April:","May:","June:","July","August","September","October","November","December")
for (i in seq(4,12)){
  print(names[i-3])
  jan <- seq(i-3,311,12)
  jan_fit <- fgev(temp_beznau$max[jan], 
                  method = "Nelder-Mead")
  par(mfrow=c(2,2))
  plot(jan_fit)
  print("Return level 100:")
  print(return_level(jan_fit,100))
  print("Return level 1000:")
  print(return_level(jan_fit,1000))
  print("Return level 10000:")
  print(return_level(jan_fit,10000))
}

# We now examine the asymptotic dependance by using chi-chibar plots
par(mfrow=c(1,2))
label_chi = sprintf("Chi plot")
label_chibar = sprintf("Chi Bar plot")
chiplot(temp, main1 = label_chi, which=1)
abline(a=0,b=0,col="red")
chiplot(temp, main2 = label_chibar, which=2)
abline(a=0,b=0,col="red")
