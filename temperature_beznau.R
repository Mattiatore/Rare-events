# loading data
library(ismev);
library(evd);
load(file="Project2020.Rdata")

# allow to print entire dataset
options(max.print=2000)

# clear seasonality, maxima are during summer
plot(ts(temp_beznau$max))

#remove 2013 since only two months?
temp_beznau=temp_beznau[-2:-1,]

# extract annual maxima
annual<-integer((nrow(temp_beznau)-11)/12+1)
annual[1]<-max(temp_beznau$max[1:11])
annual[2:length(annual)]<-sapply(split(temp_beznau$max[12:nrow(temp_beznau)], rep(1:(length(temp_beznau$max[12:nrow(temp_beznau)])/12), each=12)), max)

# fitting GEV for stationary case
monthly_fit = gev.fit(temp_beznau$max)
annual_fit = gev.fit(annual)
K <- 0
n <- length(annual)
AIC_annual <- 2*annual_fit$nllh+2*(2*K+3) +2*(2*K+3)*(2*K+4)/(n-2*K-4)
n <- nrow(temp_beznau)
AIC_mon <- 2*monthly_fit$nllh+2*(2*K+3) +2*(2*K+3)*(2*K+4)/(n-2*K-4)

# model with trig functions
# K from 1 to 9
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
    t[,i+1] <- sin(2*(i+1)*pi*diff/365.25)
  }
  # cos
  for (i in seq(K)) {
    t[,i+K+1] <- cos(2*(i+1)*pi*diff/365.25)
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

# plot for best K=6
K <- 6
t <- matrix(ncol=K*2+1,nrow=nrow(temp_beznau))
t[,1] <- diff/(100*365.25)
# sine 
for (i in seq(K)) {
  t[,i+1] <- sin(2*(i+1)*pi*diff/365.25)
}
# cos
for (i in seq(K)) {
  t[,i+K+1] <- cos(2*(i+1)*pi*diff/365.25)
}
monthly_trig <- gev.fit(temp_beznau$max,ydat=t,mul=c(1:(2*K+1)))

# plot diagnostic for constant
gev.diag(monthly_trig)

if((monthly_trig$mle[1]-1.96*monthly_trig$se[1])<0){
  print("There is no clear trend because u_1 is not different from 0 at the 95% CI")
} else{
  print("There is a trend because u_1 is different from 0 at the 95% CI")
}

# there is seasonal variation because some coefficients of the sin and cos are different from 0 at the 95% CI
