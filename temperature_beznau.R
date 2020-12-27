# loading data
library(ismev);
library(evd);
load(file="Project2020.Rdata")

# allow to print entire dataset
options(max.print=2000)

# clear seasonality, maximum during summer
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

# linear model GEV for monthly maxima
t=matrix(ncol=1,nrow=nrow(temp_beznau))
t[,1]=seq(1,nrow(temp_beznau),1)
gev.fit(temp_beznau$max,ydat=t,mul=1)

# model with trig functions?
t=matrix(ncol=9,nrow=nrow(temp_beznau))
# choose t0 as 155
t[,1] <- (seq(1,nrow(temp_beznau),1)-155)/(100+365.25)
# sine 
for (i in seq(1,4,1)) {
  t[,i+1] <- sin(2*(i+1)*pi*(seq(1,nrow(temp_beznau),1)-155)/365.25)
}
# cos
for (i in seq(1,4,1)) {
  t[,i+5] <- sin(2*(i+5)*pi*(seq(1,nrow(temp_beznau),1)-155)/365.25)
}
gev.fit(temp_beznau$max,ydat=t,mul=c(1:9))

