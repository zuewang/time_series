#load libraries
library("forecast")
library("tseries")

set.seed(42)

#load data from .csv file

ErieLevel <- read.csv('monthly-lake-erie-levels-1921-19.csv',stringsAsFactors = FALSE)
names(ErieLevel)
tail(ErieLevel)
ErieLevel <- read.csv('monthly-lake-erie-levels-1921-19.csv',header = TRUE, sep = ",")
#get month and room nights columns
names(ErieLevel)[2]
month <- ErieLevel$"Month"
level <- ErieLevel$"Monthly.Lake.Erie.Levels.1921...1970."

level_month <- ts(level, start=c(1921,1), end=c(1970,12), frequency=12)
plot(level_month, main="Raw data",xlab="Time",ylab="Monthly Lake Erie Levels")

acf(level_month, main="ACF of raw data")
#log?
sp_raw <- spectrum(level_month,main="Periodogram of raw data")
sp_raw
maxx = max(sp_raw)
abline(v=1/24, col="green",lty = 2)

#KPSS Level = 1.7146, Truncation lag parameter = 5, p-value < 0.01
kpss.test(level_month)

stl_raw <- stl(level_month,  s.window="periodic", robust=TRUE) #t.window=15,?
seasonal <- stl_raw$time.series[,1]
trend <- stl_raw$time.series[,2]
remainder <- stl_raw$time.series[,3]
kpss.test(remainder)
kpss.test(trend)
plot(seasonal)
plot(stl_raw)



#First difference plot, acf and pacf
first_diff = diff(level_month)
plot(first_diff, main="First difference",xlab="Time",ylab="Difference of Room nights by month")
spectrum(first_diff,main="Periodogram of first difference")
acf(first_diff, main="First difference")
pacf(first_diff, main="First difference")

#KPSS Level = 0.0093312, Truncation lag parameter = 5, p-value > 0.1
kpss.test(first_diff)


fit_arma <- arma(level_month, order=c(2,0))
plot(fit_arma)
summary(fit_arma)

#SARIMA
fit1 <- arima(level_month,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12))
tsdiag(fit2,font.main=3)

fit2 <- auto.arima(level_month)
tsdiag(fit2,font.main=3)
cpgram(fit2$resid,main="Cummulative Periodogram")
qqnorm(level_month)
qqline(level_month, col = 2)
qqnorm(fit2$resid)
qqline(fit2$resid, col = 2)
summary(fit2)
AIC(fit2)
BIC(fit2)
plot(forecast(fit2,h=20))
plot(fitted(fit2))





# #monthly room nights
# night_monthly <- ts(night,start=c(1921,1),frequency=12)
# 
# #Raw data plot, acf and periodogram spectrum
# plot(night_monthly, main="Raw data",xlab="Time",ylab="Room nights")
# acf(night_monthly, main="ACF of raw data")
# spectrum(night_monthly,main="Periodogram of raw data")
# 
# #First difference plot, acf and pacf
# first_diff = diff(night_monthly)
# plot(first_diff, main="First difference",xlab="Time",ylab="Difference of Room nights by month")
# acf(first_diff, main="First difference")
# pacf(first_diff, main="First difference")


