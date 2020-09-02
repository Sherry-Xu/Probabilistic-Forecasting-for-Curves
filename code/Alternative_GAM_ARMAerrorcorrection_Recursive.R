rm(list = ls())
library(magrittr)
library(yarrr)   ###for colors
library(mgcv)
library(foreach)
library(doParallel)

################################################################################
############# Train GAM recursively for all the out-of-sample data
################################################################################
## Read data
Data <- readRDS(file="../data/Data_RTE.RDS") #
Data$Month = strftime(Data$DateD,"%m") 
Data$Day =  strftime(Data$DateD,"%d") 
Data[Data$WeekDays == "Monday",'DayofWeek'] = 'Mon'
Data[Data$WeekDays %in% c("Tuesday","Wednesday","Thursday"),'DayofWeek'] = 'Tue-Thu'
Data[Data$WeekDays == "Friday",'DayofWeek'] = 'Fri'
Data[Data$WeekDays == "Saturday",'DayofWeek'] = 'Sat'
Data[Data$WeekDays == "Sunday",'DayofWeek'] = 'Sun'
summary(Data)

# Predict at noon
lag = 12*2
Data[1:(dim(Data)[1]-lag), c("Load","Forecast_RTE","Load.48","Load.336","Temp")] =
  Data[(lag+1):dim(Data)[1],c("Load","Forecast_RTE","Load.48","Load.336","Temp")]
Data = Data[1:(dim(Data)[1]-48),]
print(dim(Data)[1]/48)

# Remove holidays
sel <- which(Data$BH==1)
sel2 <- c(sel-48, sel, sel+48)
sel3 = c(sel2, which(Data$Summer_break!=0), which(Data$Christmas_break!=0))
Inval <- pmax(c(sel2, which(Data$Summer_break!=0), which(Data$Christmas_break!=0)),1)

# Remove Christmas_break
length(Inval)/48
SummerBreak = Data[(Data$Summer_break != 0) & (Data$BH!=1),]
SummerBreak$Month = "SummerBreak(08)"
ChristmasBreak = Data[Data$Christmas_break != 0 & (Data$BH!=1),]
ChristmasBreak$Month = "ChristmasBreak"
Data <- Data[-Inval,]
Data[(Data$Month == '08'),'Month']  = '07'
Data = rbind(Data, SummerBreak)  #ChristmasBreak
Data = Data[order(Data$Date),]

DataSelect = Data 
date_fin <- as.POSIXct(strptime("2018-12-31 23:30:00", "%Y-%m-%d %H:%M:%S"), tz="UTC")
date_fin2 <- as.POSIXct(strptime("2019-12-31 23:30:00", "%Y-%m-%d %H:%M:%S"), tz="UTC")
Datelist <- sort(unique(DataSelect[(DataSelect$Date>date_fin) & (DataSelect$Date<date_fin2),]$DateD))
Datelist = c(as.Date(date_fin),Datelist)
print(Datelist)

# Parallel setting
corenum = 4
myCluster <- makeCluster(corenum, # number of cores to use
                         type = "PSOCK") # type of cluster
clusterEvalQ(myCluster,library("mgcv","magrittr","yarrr"))
registerDoParallel(myCluster)

file = paste("GAM", "result.csv",sep = '_')
df_table = data.frame()
df_table_pararell = foreach(k = 1:(length(Datelist)-1),.combine = 'rbind') %dopar% {
  
  print(k)
  print(Datelist[k+1])
  
  date_fin_select <- as.POSIXct(strptime(paste(Datelist[k], " 23:30:00"), "%Y-%m-%d %H:%M:%S"), tz="UTC")
  date_fin_select2 <- as.POSIXct(strptime(paste(Datelist[k+1], " 23:30:00"), "%Y-%m-%d %H:%M:%S"), tz="UTC")
  
  Dataa <- DataSelect[DataSelect$Date<=date_fin_select,]
  Datab <- DataSelect[(DataSelect$Date>date_fin_select) & (DataSelect$Date<=date_fin_select2),]
  
  for(i in  c (0:47)){
    
    print(i)
    sel <- which(Dataa$tod==i)
    Data0 <- Dataa[sel, ]
    sel <- which(Datab$tod==i)
    Data1 <- Datab[sel, ]
    eq1_mean <- Load ~ s(as.numeric(Date), k=3)  + WeekDays:DLS  + s(toy, k = 20, bs = 'cc')+ s(as.numeric(Date),Temp, k=3) +
      s(Temp_s95, k=5) + s(Temp_s99, k=5) + s(Temp_s99_min, Temp_s99_max) + Load.48:WeekDays + Load.336 
    g1 <- gam(eq1_mean, data=Data0, family="gaussian")
    
    forecast <- predict(g1, newdata=Data1)
    df = data.frame(list("Date" = as.Date(Data1$Date), "Forecast" = forecast, "Load" = Data1$Load,"tod"=Data1$tod))
    
    forecast_train <- predict(g1, newdata=Data0)
    df_train = data.frame(list("Date" = as.Date(Data0$Date), "Forecast" = forecast_train, "Load" = Data0$Load,"tod"=Data0$tod ))
    
    residuals =  df_train$Load - df_train$Forecast
    
    resi_mean = mean(residuals)
    cat("residuals",resi_mean)
    resi_std = sd(residuals)
    
    lq = qnorm(0.05,resi_mean,resi_std)
    uq = qnorm(0.95,resi_mean,resi_std)
    
    df$lowerquantile = df$Forecast + rep(lq,length(df$Forecast))
    df$upperquantile = df$Forecast + rep(uq,length(df$Forecast))
    df_table = rbind(df_table,df)
    write.table(df,file,sep=',',col.names=FALSE,append = TRUE)
    df
  }
}

colnames = c("index","Date", "Forecast" , "Load","tod", "lowerquantile","upperquantile")

df_table = read.csv(file, header = FALSE)
colnames(df_table) = colnames
df_table = data.frame(df_table)
head(df_table)
# N = dim(df_table_train)[1]/48
Ns = dim(df_table)[1]/48
cat(Ns)

## Obtain the results without error correction
library(dplyr)
df_table = arrange(df_table, Date, tod)

MAE = mean(abs(df_table$Forecast - df_table$Load))
RMSE  = sqrt(mean((df_table$Forecast - df_table$Load)**2))
MAPE = mean(abs((df_table$Forecast - df_table$Load)/df_table$Load))
cat(MAE,RMSE,MAPE)

df_table$WithinRange = (df_table$Load >= df_table$lowerquantile) & (df_table$Load <= df_table$upperquantile)
CurveWithinRange = matrix(df_table$WithinRange, ncol = 48,byrow = TRUE)
df_table$CurveWithinRange = rep(rowSums(CurveWithinRange) == 48,each = 48)

PointCR = mean(df_table$WithinRange)
AvL = mean(df_table$upperquantile - df_table$lowerquantile)
CurveCR = mean(rowSums(CurveWithinRange) == 48)
cat("CurveCR:",CurveCR, "AvL:", AvL, "PointCR:",PointCR)

df_table$month =  format(as.Date(df_table$Date),"%m")

df_table %>%
  summarise("MAPE" = mean(abs((Forecast-Load)/Load)),
            "CurveCR" = mean(CurveWithinRange),
            "PointCR" = mean(WithinRange),
            "AvLChi2" = mean(upperquantile-lowerquantile))

df_table %>%
  group_by(month) %>%
  summarise("MAPE" = mean(abs((Forecast-Load)/Load)),
            "CurveCR" = mean(CurveWithinRange),
            "PointCR" = mean(WithinRange),
            "AvLChi2" = mean(upperquantile-lowerquantile))


################################################################################
############# Compute in sample residuals with block CV
################################################################################
library(ranger)
library(magrittr)
library(yarrr)
library(mgcv)
library(opera)
library(forecast)
source("cvblock.R")
source('buildBlock.R')

# Read data
Data <- readRDS(file="../data/Data_RTE.RDS") #
Data$Month = strftime(Data$DateD,"%m") 
Data$Day =  strftime(Data$DateD,"%d") 
Data[Data$WeekDays == "Monday",'DayofWeek'] = 'Mon'
Data[Data$WeekDays %in% c("Tuesday","Wednesday","Thursday"),'DayofWeek'] = 'Tue-Thu'
Data[Data$WeekDays == "Friday",'DayofWeek'] = 'Fri'
Data[Data$WeekDays == "Saturday",'DayofWeek'] = 'Sat'
Data[Data$WeekDays == "Sunday",'DayofWeek'] = 'Sun'
summary(Data)

# Predict at noon
lag = 12*2
Data[1:(dim(Data)[1]-lag), c("Load","Forecast_RTE","Load.48","Load.336","Temp")] =
  Data[(lag+1):dim(Data)[1],c("Load","Forecast_RTE","Load.48","Load.336","Temp")]
Data = Data[1:(dim(Data)[1]-48),]
print(dim(Data)[1]/48)

# Remove holidays
sel <- which(Data$BH==1)
sel2 <- c(sel-48, sel, sel+48)
sel3 = c(sel2, which(Data$Summer_break!=0), which(Data$Christmas_break!=0))
Inval <- pmax(c(sel2, which(Data$Summer_break!=0), which(Data$Christmas_break!=0)),1)

# Remove Christmas_break
length(Inval)/48
SummerBreak = Data[(Data$Summer_break != 0) & (Data$BH!=1),]
SummerBreak$Month = "SummerBreak(08)"
ChristmasBreak = Data[Data$Christmas_break != 0 & (Data$BH!=1),]
ChristmasBreak$Month = "ChristmasBreak"
Data <- Data[-Inval,]
Data[(Data$Month == '08'),'Month']  = '07'
Data = rbind(Data, SummerBreak)  #ChristmasBreak
Data = Data[order(Data$Date),]

date_fin <- as.POSIXct(strptime("2018-12-31 23:30:00", "%Y-%m-%d %H:%M:%S"), tz="UTC")
date_fin2 <- as.POSIXct(strptime("2019-12-31 23:30:00", "%Y-%m-%d %H:%M:%S"), tz="UTC")
Dataa <- Data[Data$Date<=date_fin,] 
Datab <- Data[(Data$Date>date_fin)&(Data$Date<=date_fin2),] 

g1 <- list()

for(i in  c (0:47))
{
  sel <- which(Dataa$tod==i)
  Data0 <- Dataa[sel, ]
  sel <- which(Datab$tod==i)
  Data1 <- Datab[sel, ]
  
  eq1_mean <- Load ~ s(as.numeric(Date), k=3)  + WeekDays:DLS  + s(toy, k = 20, bs = 'cc')+s(as.numeric(Date),Temp, k=3) +
    s(Temp_s95, k=5) + s(Temp_s99, k=5) + s(Temp_s99_min, Temp_s99_max) + Load.48:WeekDays +Load.336 
  g1[[i+1]] <- gam(eq1_mean , data=Data0, family="gaussian")
  #g1[[i+1]]$forecast <- predict(g1[[i+1]], newdata=Data1)
  
  bl <- buildBlock(Nblock=10, data=Data0)
  forecastCV <- lapply(bl, forecastBlock, formula=eq1_mean, data=Data0)
  forecastCV <- unlist(forecastCV)
  g1[[i+1]]$residualsCV <- Data0$Load-forecastCV
  print(i)
}

res <- list()
res$g1 <- g1
saveRDS(res, file="gam_with_blockCV.RDS")

################################################################################
############# ARMA Error correction
################################################################################
#### get the results of the GAM with block CV residuals
res <- readRDS(file="gam_with_blockCV.RDS")

#### compute the CV block error on the train set, for the estimation of the ARIMA intraday correction
resCV <- array(0, dim=nrow(Dataa))
for(i in c(1:48))
{
  sel <- which(Dataa$tod==i-1)
  resCV[sel] <- res$g1[[i]]$residualsCV
}
#### get the results of the GAM recursive forecasting
file = paste("GAM", "result.csv",sep = '_')
colnames = c("index","Date", "Forecast" , "Load","tod", "lowerquantile","upperquantile")
df_table = read.csv(file, header = FALSE)
colnames(df_table) = colnames
df_table = data.frame(df_table)
df_table = arrange(df_table, Date, tod)
head(df_table)
df_table$Date = as.Date(df_table$Date)
dim(df_table)[1]/48

#res.test <- Datab$Load-gam.prev
res.test <- df_table$Load - df_table$Forecast

####"box jenkins" like analysis
plot(resCV, type='l')
acf(resCV)
pacf(resCV)
par(mfrow=c(1,2))
acf(diff(resCV))
pacf(diff(resCV))

plot(diff(resCV), type='l')

sqrt(mean(resCV^2))
sqrt(mean(res.test^2))

####fit the ARIMA model
resCV.ts <- ts(resCV)
ar <- Arima(resCV.ts, order=c(4,0,0), method = c("CSS"))
ar

#### do the rolling short term correction on the test set
ar.forecast <- forecast(ar, h=48, level=c(90))
pred.ar <- ar.forecast$mean
pred.ar.lower <- ar.forecast$lower ###upper bound of the 90% quantile interval
pred.ar.upper <- ar.forecast$upper

for(i in c(1:(nrow(df_table)/48-1))){
  res.pred <- c(resCV.ts, res.test[((i-1)*48+1):(i*48)])
  res.pred <- ts(res.pred)
  refit <- Arima(res.pred, model=ar)
  ar.forecast <- forecast(refit, h=48, level=c(90))
  pred.ar <- c(pred.ar, ar.forecast$mean)
  pred.ar.lower <- c(pred.ar.lower, ar.forecast$lower)
  pred.ar.upper <- c(pred.ar.upper, ar.forecast$upper)
}

gam.prev = df_table$Forecast
gam.intraday <- gam.prev+pred.ar


####check
a <- 48*7*2
b <- a+48

par(mfrow=c(1,1))
plot(Datab$Load[a:b], type='l')
lines(gam.prev[a:b], col='red')
lines(gam.intraday[a:b], col='blue')

# Get the results with error correction
lowerquantile = gam.prev+pred.ar.lower
upperquantile = gam.prev+pred.ar.upper
df_table = data.frame(list("Date" = as.Date(df_table$Date), "Forecast_gam" = gam.prev,
                           "Forecast" = gam.intraday,
                           "Load" = df_table$Load,"tod"=df_table$tod,
                           "lowerquantile" = lowerquantile,
                           "upperquantile" = upperquantile))

head(df_table)
dim(df_table)[1]/48

library(dplyr)
df_table = arrange(df_table, Date, tod)

MAE = mean(abs(df_table$Forecast - df_table$Load))
RMSE  = sqrt(mean((df_table$Forecast - df_table$Load)**2))
MAPE = mean(abs((df_table$Forecast - df_table$Load)/df_table$Load))
cat(MAE,RMSE,MAPE)

df_table$WithinRange = (df_table$Load >= df_table$lowerquantile) & (df_table$Load <= df_table$upperquantile)
CurveWithinRange = matrix(df_table$WithinRange, ncol = 48,byrow = TRUE)
df_table$CurveWithinRange = rep(rowSums(CurveWithinRange) == 48,each = 48)

PointCR = mean(df_table$WithinRange)
AvL = mean(df_table$upperquantile - df_table$lowerquantile)
CurveCR = mean(rowSums(CurveWithinRange) == 48)
cat("CurveCR:",CurveCR, "AvL:", AvL, "PointCR:",PointCR)


df_table$month =  format(as.Date(df_table$Date),"%m")

df_table %>%
  summarise("MAPE" = mean(abs((Forecast-Load)/Load)),
            "CurveCR" = mean(CurveWithinRange),
            "PointCR" = mean(WithinRange),
            "AvLChi2" = mean(upperquantile-lowerquantile))

df_table %>%
  group_by(month) %>%
  summarise("MAPE" = mean(abs((Forecast-Load)/Load)),
            "CurveCR" = mean(CurveWithinRange),
            "PointCR" = mean(WithinRange),
            "AvLChi2" = mean(upperquantile-lowerquantile))

(0.0124+0.00950+0.0106+0.0145)/4
