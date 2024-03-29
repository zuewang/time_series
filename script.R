#load libraries
library("forecast")
library("tseries")
set.seed(42)
par(mar=c(5.1,4.1,4.1,2.1))
Hotel <- read.csv('monthly-data-relating-to-hotels-.csv',header = TRUE, sep = ",", stringsAsFactors = FALSE)
#get month, room nights, takings columns
len <- length(Hotel$"Month")
#the last one is invalid
month <- Hotel$"Month"[1:len-1]
nights <- Hotel$"Room_nights"[1:len-1]
takings <- Hotel$"Takings"[1:len-1]

nights <- as.integer(nights)
nights_month <- ts(nights, start=c(1980,1), end=c(1995,6), frequency=12)
takings_month <- ts(takings, start=c(1980,1), end=c(1995,6), frequency=12)

#plot raw data
par(mfrow=c(2,1))
plot(nights_month, xlab="Time", ylab="Number", main="Monthly Room Nights")
plot(takings_month, xlab="Time", ylab="Number", main="Monthly Takings")

#plot takings vs nights 
par(mfrow=c(1,1))
cor(nights_month, takings_month)
cor.test(nights_month, takings_month)
plot(nights_month,takings_month, main="Monthly takings vs Monthly Room_nights")

# ACF PACF
par(mfrow=c(1,2))
#acf(nights_month, main="ACF of Room_nights")
#pacf(nights_month, main="PACF of Room_nights")
spectrum(nights_month,main="Periodogram of Room_nights")
#acf(takings_month, main="ACF of Takings")
#pacf(takings_month, main="PACF of Takings")
spectrum(takings_month,main="Periodogram of Takings")

#KPSS test for stationarity of raw data
kpss.test(nights_month)
kpss.test(takings_month)

#STL decomposition
par(mfrow=c(1,1))
stl_nights <- stl(nights_month,  s.window="periodic")#, robust=TRUE) t.window=15,?
plot(stl_nights)
stl_takings <- stl(takings_month,  s.window="periodic")
plot(stl_takings)

#First difference plot, acf and pacf
dy_n = diff(nights_month)
dy_t = diff(takings_month)
par(mfrow=c(2,4))
plot(dy_n, main="First difference of Room_nights (dy_n)",xlab="Time",ylab="First difference")
acf(dy_n, lag.max = 48, main="ACF of dy_n")
pacf(dy_n, lag.max = 48, main="PACF of dy_n")
qqnorm(dy_n/sd(dy_n), main="Q-Q plot of standardized dy_n"); abline(a=0,b=1)
plot(dy_t, main="First difference of Takings (dy_t)",xlab="Time",ylab="First difference")
acf(dy_t, lag.max = 48, main="ACF of dy_t")
pacf(dy_t, lag.max = 48, main="PACF of dy_t")
qqnorm(dy_t/sd(dy_t), main="Q-Q plot of standardized dy_t"); abline(a=0,b=1)

#KPSS test for stationarity of first differences
kpss.test(dy_n)
kpss.test(dy_t)

#fit SARIMA models:
#(0,0,1)(0,1,1), (1,0,1)(0,1,1), (0,0,1)(0,1,2), (1,0,1)(0,1,2)
#for room_nights
fit_dyn1 <- arima(dy_n,order=c(0,0,1),seasonal=list(order=c(0,1,1),period=12))
fit_dyn2 <- arima(dy_n,order=c(1,0,1),seasonal=list(order=c(0,1,1),period=12))
fit_dyn3 <- arima(dy_n,order=c(0,0,1),seasonal=list(order=c(0,1,2),period=12))
fit_dyn4 <- arima(dy_n,order=c(1,0,1),seasonal=list(order=c(0,1,2),period=12))
par(mfrow=c(2,2))
cpgram(fit_dyn1$resid, main='Cummulative periodogram for ARIMA(0,0,1)(0,1,1)[12]')
cpgram(fit_dyn2$resid, main='Cummulative periodogram for ARIMA(1,0,1)(0,1,1)[12]')
cpgram(fit_dyn3$resid, main='Cummulative periodogram for ARIMA(0,0,1)(0,1,2)[12]')
cpgram(fit_dyn4$resid, main='Cummulative periodogram for ARIMA(1,0,1)(0,1,2)[12]')
qqnorm(fit_dyn1$resid/sd(fit_dyn1$resid), main='Normal Q-Q Plot for standard residuals of ARIMA(0,0,1)(0,1,1)[12]'); abline(a = 0, b = 1)
qqnorm(fit_dyn2$resid/sd(fit_dyn2$resid), main='Normal Q-Q Plot for standard residuals of ARIMA(1,0,1)(0,1,2)[12]'); abline(a = 0, b = 1)
qqnorm(fit_dyn3$resid/sd(fit_dyn3$resid), main='Normal Q-Q Plot for standard residuals of ARIMA(0,0,1)(0,1,1)[12]'); abline(a = 0, b = 1)
qqnorm(fit_dyn4$resid/sd(fit_dyn4$resid), main='Normal Q-Q Plot for standard residuals of ARIMA(1,0,1)(0,1,2)[12]'); abline(a = 0, b = 1)
tsdiag(fit_dyn1)
tsdiag(fit_dyn2)
tsdiag(fit_dyn3)
tsdiag(fit_dyn4)
summary(fit_dyn1)
summary(fit_dyn2)
summary(fit_dyn3)
summary(fit_dyn4)
summary(fit_dyn1)#,fit_dyn2,fit_dyn3,fit_dyn4)
AIC(fit_dyn1,fit_dyn2,fit_dyn3,fit_dyn4)
BIC(fit_dyn1,fit_dyn2,fit_dyn3,fit_dyn4)
#for takings
fit_dyt1 <- arima(dy_t,order=c(0,0,1),seasonal=list(order=c(0,1,1),period=12))
fit_dyt2 <- arima(dy_t,order=c(1,0,1),seasonal=list(order=c(0,1,1),period=12))
fit_dyt3 <- arima(dy_t,order=c(0,0,1),seasonal=list(order=c(0,1,2),period=12))
fit_dyt4 <- arima(dy_t,order=c(1,0,1),seasonal=list(order=c(0,1,2),period=12))
par(mfrow=c(2,2))
cpgram(fit_dyt1$resid, main='Cummulative periodogram for ARIMA(0,0,1)(0,1,1)[12]')
cpgram(fit_dyt2$resid, main='Cummulative periodogram for ARIMA(1,0,1)(0,1,1)[12]')
cpgram(fit_dyt3$resid, main='Cummulative periodogram for ARIMA(0,0,1)(0,1,2)[12]')
cpgram(fit_dyt4$resid, main='Cummulative periodogram for ARIMA(1,0,1)(0,1,2)[12]')
qqnorm(fit_dyt1$resid/sd(fit_dyt1$resid), main='Normal Q-Q Plot for standard residuals of ARIMA(0,0,1)(0,1,1)[12]'); abline(a = 0, b = 1)
qqnorm(fit_dyt2$resid/sd(fit_dyt2$resid), main='Normal Q-Q Plot for standard residuals of ARIMA(1,0,1)(0,1,1)[12]'); abline(a = 0, b = 1)
qqnorm(fit_dyt3$resid/sd(fit_dyt3$resid), main='Normal Q-Q Plot for standard residuals of ARIMA(0,0,1)(0,1,2)[12]'); abline(a = 0, b = 1)
qqnorm(fit_dyt4$resid/sd(fit_dyt4$resid), main='Normal Q-Q Plot for standard residuals of ARIMA(1,0,1)(0,1,2)[12]'); abline(a = 0, b = 1)
tsdiag(fit_dyt1)
tsdiag(fit_dyt2)
tsdiag(fit_dyt3)
tsdiag(fit_dyt4)
summary(fit_dyt1)
summary(fit_dyt2)
summary(fit_dyt3)
summary(fit_dyt4)
AIC(fit_dyt1,fit_dyt2,fit_dyt3,fit_dyt4)
BIC(fit_dyt1,fit_dyt2,fit_dyt3,fit_dyt4)

# one-step prediction
train_len <- as.integer((len-1)/2)
predict_len <- len-1-train_len
nights_train <- ts(nights, start=c(1980,1), end=c(1987,9), frequency=12)
takings_train <- ts(takings, start=c(1980,1), end=c(1987,9), frequency=12)
fit_n_train <- arima(nights_train,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12))
fit_t_train <- arima(takings_train,order=c(0,1,1),seasonal=list(order=c(0,1,2),period=12))
fore_n <- forecast(fit_n_train, predict_len)
fore_t <- forecast(fit_t_train, predict_len)

par(mfrow=c(1,1))
# Plot Room_nights
plot(nights_month,main='Forecasts from ARIMA(0,1,1)(0,1,1)[12]',ylab='Monthly Room_nights')
lines(fore_n$fitted, col='red')
lines(fore_n$mean, col='blue')
# 95% confidence interval
lines(fore_n$upper[,2], lty=2, col='green')
lines(fore_n$lower[,2], lty=4, col='green')
legend('topleft',
       c('Raw data','Fitted', 'Predicted mean', '95% CI upper bound', '95% lower bound'), 
       lty=c(1,1,1,2,4), col=c('black','red','blue','green','green'))
# Plot Takings
plot(takings_month,main='Forecasts from ARIMA(0,1,1)(0,1,2)[12]',ylab='Monthly Takings')
lines(fore_t$fitted, col='red')
lines(fore_t$mean, col='blue')
# 95% confidence interval
lines(fore_t$upper[,2], lty=2, col='green')
lines(fore_t$lower[,2], lty=4, col='green')
legend('topleft',
       c('Raw data','Fitted', 'Predicted mean', '95% CI upper bound', '95% lower bound'), 
       lty=c(1,1,1,2,4), col=c('black','red','blue','green','green'))

# Plot forcasts with total data fitted
plot(forecast(arima(nights_month,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12))), ylab='Monthly Room_nights')
plot(forecast(arima(takings_month,order=c(0,1,1),seasonal=list(order=c(0,1,2),period=12))), ylab='Monthly Takings')






#rectify Error in plot.new() : figure margins too large
#par(mar=c(1,1,1,1))
#par(mar = rep(2, 4))

#load data from .csv file
#Hotel <- read.csv('monthly-data-relating-to-hotels-.csv',stringsAsFactors = FALSE)
names(Hotel)
tail(Hotel)



# uncertainty 


# q-q plot


#0  0  1  0  1  1 12
fit1 <- auto.arima(dy_n)
arimaorder(fit1)
AIC(fit1)
#1  0  1  0  1  1 12 not good as the above one
fit11 <- auto.arima(nights_month)
fit11 <- arima(nights_month,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12))
arimaorder(fit11)
AIC(fit11)
# 0  0  1  0  1  2 12
fit2 <- auto.arima(dy_t)
arimaorder(fit2)
AIC(fit2)

#0  1  1  0  1  2 12
fit22 <- auto.arima(takings_month)
arimaorder(fit22)
AIC(fit22)

fit_dy1_1 <- arima(first_diff,order=c(0,0,0),seasonal=list(order=c(2,0,0),period=12), transform.pars = TRUE)
arimaorder(fit_dy1_1)
#fit_dy1_sarima <- sarima(first_diff, p = 2, d = 0, q = 1, details = TRUE)
plot(fit_dy1_1)
AIC(fit_dy1_1)
summary(fit_dy1_1)
tsdiag(fit_dy1_1,font.main=3)

fit_arma <- arma(nights_month, order=c(2,0))
plot(fit_arma)
summary(fit_arma)


#SARIMA
fit1 <- arima(nights_month,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12))
tsdiag(fit1,font.main=3)

fit2 <- auto.arima(nights_month)
arimaorder(fit2)
tsdiag(fit2,font.main=3)
cpgram(fit2$resid,main="Cummulative Periodogram")
qqnorm(nights_month)
qqline(nights_month, col = 2)
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




#log?
sp_raw <- spectrum(nights_month,main="Periodogram of raw data")
sp_raw
maxx = max(sp_raw)
abline(v=1/24, col="green",lty = 2)

seasonal <- stl_nights$time.series[,1]
trend <- stl_nights$time.series[,2]
stl_res_nights <- stl_nights$time.series[,3]
kpss.test(stl_res_nights)
fit_stl_res <- auto.arima(stl_res)
arimaorder(fit_stl_res)
kpss.test(trend)
tsdiag(fit_stl_res)

plot(seasonal)


#regression
ts1_a <- lm(ErieLevel ~ month + fourier(ErieLevel, K = 1/12))
forecast::Acf(resid(ts1_a), main = "Correlogram of residuals")


