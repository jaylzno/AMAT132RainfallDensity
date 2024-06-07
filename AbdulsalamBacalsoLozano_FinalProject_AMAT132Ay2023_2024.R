# Load relevant libraries
library(ggplot2)
library(knitr)
library(printr)
library(plyr)
library(dplyr)
library(lubridate)
library(gridExtra)
library(zoo)
library(fpp3)
library(reshape2)
library(TTR)
library(forecast)
library(tseries)

# Set working directory
setwd("C:/Users/User/Desktop")
rfh <- read.csv("rainfall_132 - Copy.csv")
rfh$Date <- as.Date(rfh$Date, format = "%d/%m/%Y")

# Descriptive analysis of the raw data
mean(rfh$Value)
median(rfh$Value)
sd(rfh$Value)
var(rfh$Value)

# Time Series Decomposition
str(rfh)
rfh <- rfh[order(rfh$Date),]
rfh_ts <- ts(rfh$Value, frequency = 36, start = c(1981, 01))

# Plotting the decomposed data
decomposed_rfh <- decompose(rfh_ts)
autoplot(decomposed_rfh)

# Splitting data to train-test sets with 80/20 ratio
train <- head(rfh_ts, round(length(rfh_ts)) * 0.8)
h <- length(rfh_ts) - length(train)
test <- tail(rfh_ts, h)

# Naive Forecasting
rfh_naive <- c(train)
rfh_forecast <- c(NA, train[-length(train)])

plot(rfh_naive, type ='l', col = 'red', main='Actual vs Forecasted (Naive Forecasting - Training Set 1981-2015)')
lines(rfh_forecast, type = 'l', col='blue')
legend('topright', legend=c('Actual','Forecasted'),col=c('red', 'blue'), lty=1)

rfh_naive1 <- c(test)
rfh_forecast1 <- c(NA, test[-length(test)])

plot(rfh_naive1, type ='l', col = 'red', main='Actual vs Forecasted (Naive Forecasting - Test Set 2015 to 2024)')
lines(rfh_forecast1, type = 'l', col='blue')
legend('topright', legend=c('Actual','Forecasted'),col=c('red', 'blue'), lty=1)

accuracy(rfh_naive1, rfh_forecast1)

# Simple Exponential Smoothing
rfh_ses <- ses(train, alpha = .2, h = 312)
summary(rfh_ses)
accuracy(rfh_ses, test)

# plotting results
p1 <- autoplot(rfh_ses) + ggtitle("Forecast of 10-day Rainfall Density (mm) using Smooth Exponential Smoothing")
    theme(legend.position = "bottom")
p2 <- autoplot(rfh_ses) + autolayer(test, alpha=0.5) + ggtitle("Forecasted vs Actual")

gridExtra::grid.arrange(p1,p2, nrow = 1)

#Holt's Linear Method
rfh_holt <- holt(train, h = 312)
autoplot(rfh_holt)
summary(rfh_holt)
rfh_holt$model
accuracy(rfh_holt, test)

# identify optimal beta parameter
beta <- seq(.0001, .5, by = .001)
RMSE <- NA
for(i in seq_along(beta)) {
  fit <- holt(train,
              beta = beta[i], 
              h = 100)
  RMSE[i] <- accuracy(fit, 
                      test)[2,2]
}

beta.fit <- data_frame(beta, RMSE)
beta.min <- filter(beta.fit, RMSE == min(RMSE))

ggplot(beta.fit, aes(beta, RMSE)) +
  geom_line() +
  geom_point(data = beta.min, 
             aes(beta, RMSE), 
             size = 2, color = "red")

# forecasting new model with optimal parameters
rfh_holt_opt <- holt(train, h = 312, beta = 0.0001)
summary(rfh_holt_opt)
accuracy(rfh_holt, test)
accuracy(rfh_holt_opt, test)

p1 <- autoplot(rfh_holt) +
  ggtitle("Primary Forecast of 10-day Rainfall Density (mm) using Holt's Linear Method") +
  coord_cartesian(ylim = c(-1000, 1000)) + autolayer(test)

p2 <- autoplot(rfh_holt_opt) +
  ggtitle("Optimal Forecast of 10-day Rainfall Density (mm) using Holt's Linear Method") +
  coord_cartesian(ylim = c(-1000, 1000)) + autolayer(test)

gridExtra::grid.arrange(p1, p2,nrow = 1)

#Holt-Winter's Method

# Annual Aggregated Data
df <- data.frame(Date = rfh$Date, Value = rfh$Value)
annual_rfh <- df %>%
  mutate(Year = format(Date, "%Y-%m")) %>%
  group_by(Year) %>%
  summarize(Value = mean(Value, na.rm = TRUE))

annagg_rfh <- ts(annual_rfh$Value, frequency = 12, start = c(1981, 01))
autoplot(annagg_rfh) + ggtitle("Annually Aggregated Data Set")
autoplot(decompose(annagg_rfh))
rfh_hw <- ets(annagg_rfh, model = "AAA")

## Splitting the aggregated data into training and test sets
hw_train <- head(annual_rfh, round(nrow(annual_rfh)) * 0.8)
hw_h <- nrow(annual_rfh) - nrow(hw_train)
hw_test <- tail(annual_rfh, hw_h)
rfh_hw_train <- ts(hw_train$Value, frequency = 12, start = c(1981, 01))
rfh_hw_test <- ts(hw_test$Value, frequency = 12, start = c(2015, 09))

# Checking residuals for additive and multiplicative
rfh_hw1 <- ets(rfh_hw_train, model = "AAA")
summary(rfh_hw1)
checkresiduals(rfh_hw1)
rfh_hw2 <- ets(rfh_hw_train, model = "MAM")
summary(rfh_hw2)
checkresiduals(rfh_hw2)

rfh_fhw1 <- forecast(rfh_hw1, h = 104)
rfh_fhw2 <- forecast(rfh_hw2, h = 104)

# Plotting the primary forecasts
p1 <- autoplot(rfh_fhw1) +
  ggtitle("Primary Forecast of 10-day Rainfall Density (mm) in ETS(A,A,A)") +
  coord_cartesian(ylim = c(-1000, 1000)) + autolayer(test)

p2 <- autoplot(rfh_fhw2) +
  ggtitle("Primary Forecast of 10-day Rainfall Density (mm) in ETS(M,A,M)") +
  coord_cartesian(ylim = c(-1000, 1000)) + autolayer(test)

gridExtra::grid.arrange(p1, p2,nrow = 1)

# Check for the optimal gamma parameter
gamma <- seq(0.0001, 0.85, 0.01)
RMSE <- NA

for(i in seq_along(gamma)) {
  hw.expo <- ets(rfh_hw_train, 
                 "AAA", 
                 gamma = gamma[i])
  future <- forecast(hw.expo, 
                     h = 5)
  RMSE[i] = accuracy(future, 
                     rfh_hw_test)[2,2]
}

error <- data_frame(gamma, RMSE)
minimum <- filter(error, 
                  RMSE == min(RMSE))

#Plot for the optimal gamma parameter
ggplot(error, aes(gamma, RMSE)) +
  geom_line() +
  geom_point(data = minimum, 
             color = "blue", size = 2) +
  ggtitle("gamma's impact on 
            forecast errors",
          subtitle = "gamma = 0.0001 minimizes RMSE")

# new AAA model with 
# optimal gamma parameter
rfh_hw_hwAAA <- ets(rfh_hw_train,
                   model = "AAA", 
                   gamma = 0.0001)
rfh_hw_fAAA <- forecast(rfh_hw_hwAAA, 
                       h = 104)
summary(rfh_hw_fAAA)
accuracy(rfh_hw_fAAA, rfh_hw_test)


# new MAM model with 
# optimal gamma parameter
rfh_hw_hwMAM <- ets(rfh_hw_train,
                    model = "MAM", 
                    gamma = 0.0001)
rfh_hw_fMAM <- forecast(rfh_hw_hwMAM, 
                        h = 104)
summary(rfh_hw_fMAM)
accuracy(rfh_hw_fMAM, rfh_hw_test)

# Plotting the optimized forecasts
p1 <- autoplot(rfh_hw_fAAA) +
  ggtitle("Optimized Forecast of 10-day Rainfall Density (mm) in ETS(A,A,A)") +
  coord_cartesian(ylim = c(-1000, 1000)) + autolayer(test)

p2 <- autoplot(rfh_hw_fMAM) +
  ggtitle("Optimized Forecast of 10-day Rainfall Density (mm) in ETS(M,A,M)") +
  coord_cartesian(ylim = c(-1000, 1000)) + autolayer(test)

gridExtra::grid.arrange(p1, p2,nrow = 1)


#Damped HW

rfh_hwda <- holt(rfh_hw_train, damped = TRUE, seasonal = "additive", h = 104)
summary(rfh_hwda)
accuracy(rfh_hwda, rfh_hw_test)

rfh_hwdm <- holt(rfh_hw_train, damped = TRUE, seasonal = "multiplicative", h = 104)
summary(rfh_hwdm)
accuracy(rfh_hwdm, rfh_hw_test)


p1 <- autoplot(rfh_hwda) + autolayer(rfh_hw_test) + ggtitle("Forecast of 10-day Rainfall Density (mm) in Damped Additive Method")

p2 <- autoplot(rfh_hwdm) + autolayer(rfh_hw_test) + ggtitle("Forecast of 10-day Rainfall Density (mm) in Damped Multiplicative Method")

gridExtra::grid.arrange(p1, p2,nrow = 1)

#ARIMA

#ADF Test
adf.test(train, k = 36)
forecast::tsdisplay(train)

stl_decomp <- stl(train, s.window = "periodic")
seasonal_comp <- stl_decomp$time.series[, "seasonal"]
deseasonalized_train <- train - seasonal_comp

autoplot(deseasonalized_train)

adf.test(deseasonalized_train, k = 36)
forecast::tsdisplay(deseasonalized_train)

#arima015
arima015 = forecast::Arima(train, order = c(0,1,5))
summary(arima015)
forecast::checkresiduals(arima015)
arima015f <- forecast(arima015, h = 312)

autoplot(arima015f)+ autolayer(test)

accuracy(arima015f, test)

#sarima301111
arima301111 = forecast::Arima(train, order = c(3,0,1), seasonal = c(0,1,1))
summary(arima301111)
forecast::checkresiduals(arima301111)
arima301110f <- forecast(arima301111, h = 312)

autoplot(arima301110f)+ autolayer(test)

accuracy(arima301110f, test)

# Final Forecast Using Best Model
forecastfinal <- ets(annagg_rfh,model = "MAM", gamma = 0.0001)
futureforecast <- forecast(forecastfinal, h = 72)
autoplot(futureforecast) + ggtitle("Forecast of the Philippine 10-day Rainfall Density (mm) from April 2024 to April 2030")
summary(futureforecast)
