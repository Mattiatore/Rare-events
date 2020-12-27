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

# model with trig functions
# K from 1 to 9
# assign t0 as in the description 
t0 = as.Date("2000/01/01")
diff <- seq(nrow(temp_beznau))
for (i in seq(nrow(temp_beznau))) {
  str <- paste(toString(temp_beznau$year[i]),toString(temp_beznau$mon[i]),toString(temp_beznau$day[i]),sep="/")
  diff[i] <- as.Date(str)-t0
}
AIC_monthly <- seq(9)
for (K in seq(9)){
  t <- matrix(ncol=K*2+1,nrow=nrow(temp_beznau))
  t[,1] <- diff/(100+365.25)
  # sine 
  for (i in seq(K)) {
    t[,i+1] <- sin(2*(i+1)*pi*diff/365.25)
  }
  # cos
  for (i in seq(K)) {
    t[,i+K+1] <- sin(2*(i+K+1)*pi*diff/365.25)
  }
  monthly_trig <- gev.fit(temp_beznau$max,ydat=t,mul=c(1:2*K+1))
  
  # AIC = -2 * LL + 2 * p + 2*p*(p+1)/(n-p-1)
  AIC_monthly[K] <- 2*monthly_trig$nllh+2*K +2*K*(K+1)/(n-K-1)
}

