rm(list = ls())
library("forecast")
source("cvblock.R")
source('buildBlock.R')

PrepareData = function(DataSelect, feature_list, date_fin,date_fin2){
  DataTrain <- DataSelect[DataSelect$Date<=date_fin,]
  DataTest <- DataSelect[(DataSelect$Date>date_fin) & (DataSelect$Date<=date_fin2),]
  
  DataTrainY = DataTrain$Load
  DataTrainY = matrix(DataTrainY, ncol = 48, byrow = TRUE)
  DataTestY = DataTest$Load
  DataTestY = matrix(DataTestY, ncol = 48, byrow = TRUE)
  
  DataTrainX = matrix(0, nrow = dim(DataTrain)[1]/48,ncol = 0)
  DataTestX = matrix(0, nrow = dim(DataTest)[1]/48,ncol = 0)
  for (column in feature_list){
    DataTrainX1 = DataTrain[,column]
    DataTrainX1 = matrix(DataTrainX1, ncol = 48, byrow = TRUE)
    DataTrainX = cbind(DataTrainX,DataTrainX1)
    DataTestX1 = DataTest[,column]
    DataTestX1 = matrix(DataTestX1, ncol = 48, byrow = TRUE)
    DataTestX = cbind(DataTestX,DataTestX1)
  }
  return(data = list("DataTrainY" = DataTrainY, "DataTrainX" = DataTrainX, 
                     "DataTestY" = DataTestY, "DataTestX" = DataTestX))
}

Data <- readRDS(file="../newdata/Data_RTE.RDS")
Data$Month = strftime(Data$DateD,"%m") 
Data$Day =  strftime(Data$DateD,"%d") 
Data[Data$WeekDays == "Monday",'DayofWeek'] = 'Mon'
Data[Data$WeekDays %in% c("Tuesday","Wednesday","Thursday"),'DayofWeek'] = 'Tue-Thu'
Data[Data$WeekDays == "Friday",'DayofWeek'] = 'Fri'
Data[Data$WeekDays == "Saturday",'DayofWeek'] = 'Sat'
Data[Data$WeekDays == "Sunday",'DayofWeek'] = 'Sun'
summary(Data)

Data$Load.96 = 0
Data$Load.96[(48*2+1):(146064)] = Data$Load[1:(146064-48*2)]
Data = Data[(48*2+1):(dim(Data)[1]-48),]

lag = 12*2
Data[1:(dim(Data)[1]-lag), c("Load","Forecast_RTE","Load.48","Load.96","Load.336","Temp")] =
  Data[(lag+1):dim(Data)[1],c("Load","Forecast_RTE","Load.48","Load.96","Load.336","Temp")]
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
Data = rbind(Data, SummerBreak)
Data = Data[order(Data$Date),]

date_fin <- as.POSIXct(strptime("2018-12-31 23:30:00", "%Y-%m-%d %H:%M:%S"), tz="UTC")
date_fin2 <- as.POSIXct(strptime("2019-12-31 23:30:00", "%Y-%m-%d %H:%M:%S"), tz="UTC")

splitweekday = "DayofWeek"
#feature = "17T"
feature_list = c("Load.48","Load.96","Load.336","Temp") 

Data[Data$Month %in% c("06","07","09"), "Month"] = "06-09" 
Data[Data$Month %in% c("12","01","02"), "Month"] = "12-02"
Data[Data$Month %in% c("04","05"), "Month"] = "04-05"
#Data[Data$Month %in% c("03","11"), "Month"] = "03,11"

weekdays_list = unique(Data[,splitweekday])
print(weekdays_list)
months_list = sort(unique(Data$Month))
print(months_list)

month = "12-02" #ChristmasBreak"#,"02,03" #months_list[12]
weekday = weekdays_list[5]


df_table = data.frame()

#file = paste("20200701_Real", lag, feature,Alpha_list,Dhat_s,"month_result.csv",sep = '_')
for(month in months_list){
  for (weekday in weekdays_list){
    print(weekday)
    print(month)
    
    DataSelect = Data[(Data$DayofWeek == weekday) & (Data$Month == month),] #& 
    Datelist <- sort(unique(DataSelect[(DataSelect$Date>date_fin) & (DataSelect$Date<date_fin2),]$DateD))
    Datelist = c(as.Date(date_fin),Datelist)
    print(Datelist)
    
    #
    for (i in 1:(length(Datelist)-1)){
      
      print(i)
      print(Datelist[i+1])
      date_fin_select <- as.POSIXct(strptime(paste(Datelist[i], " 23:30:00"), "%Y-%m-%d %H:%M:%S"), tz="UTC")
      date_fin_select2 <- as.POSIXct(strptime(paste(Datelist[i+1], " 23:30:00"), "%Y-%m-%d %H:%M:%S"), tz="UTC")
      
      data = PrepareData(DataSelect, feature_list,date_fin_select,date_fin_select2)
      
      N = dim(data$DataTrainY)[1]
      Ns = dim(data$DataTestY)[1]
      
      cat(N,Ns,"\n")
      Y = data
      DataType = "Real"
      
      YY = Y$DataTrainY
      XX = Y$DataTrainX
      y = Y$DataTestY
      x = Y$DataTestX
      
      #if (i == 1){
      Mean = colMeans(YY)
      Meanx = colMeans(XX)
      Std = apply(YY, 2, sd)
      Stdx = apply(XX, 2, sd)
      #}
      
      YY = (YY - matrix(rep(Mean, N), nrow=N, byrow=T))
      XX = (XX - matrix(rep(Meanx, N), nrow=N, byrow=T))
      y = (y - matrix(rep(Mean, Ns), nrow=Ns, byrow=T))
      x = (x - matrix(rep(Meanx, Ns), nrow=Ns, byrow=T))

      YY = YY /matrix(rep(Std, N), nrow=N, byrow=T)
      XX = XX/matrix(rep(Stdx, N), nrow=N, byrow=T)
      y = y /matrix(rep(Std, Ns), nrow=Ns, byrow=T)
      x = x/matrix(rep(Stdx, Ns), nrow=Ns, byrow=T)
      
      YYh = matrix(0,N,48)
      yh = matrix(0,Ns,48)
      
      # SARX fitting
      for (j in 1:48){
        nrep = dim(XX)[1]
        Xm=matrix(0,N,4)
        Xm[,1] = XX[,j]
        Xm[,2] = XX[,(j+48)]
        Xm[,3] = XX[,(j+96)]
        Xm[,4] = XX[,(j+144)]
        Ym = YY[,j]
        beta <- solve(t(Xm)%*%Xm)%*%t(Xm)%*%matrix(Ym,ncol=1) 
        
        xm=matrix(0,Ns,4)
        xm[,1] = x[,j]
        xm[,2] = x[,(j+48)]
        xm[,3] = x[,(j+96)]
        xm[,4] = x[,(j+144)]
        #YYh[,j] = Xm%*%beta
        yh[,j]=xm%*%beta
      }
      
      # Compute cross-validation residuals
      for (j in 1:48){
        Xm=matrix(0,N,4)
        Xm[,1] = XX[,j]
        Xm[,2] = XX[,(j+48)]
        Xm[,3] = XX[,(j+96)]
        Xm[,4] = XX[,(j+144)]
        Ym = YY[,j]
        bl  = buildBlock(10,XX)
        for (kbl in 1:10){
        XmCV = Xm[-bl[[kbl]],]
        YmCV = Ym[-bl[[kbl]]]
        betaCV <- solve(t(XmCV)%*%XmCV)%*%t(XmCV)%*%matrix(YmCV,ncol=1) 
        YYh[bl[[kbl]],j] = Xm[bl[[kbl]],]%*%betaCV
        }
      }
      
      residuals = YY - YYh
      residuals = matrix(t(residuals),ncol=1)
      length(residuals)
      
      ####"box jenkins" like analysis
      plot(residuals, type='l')
      acf(residuals)
      pacf(residuals)
      par(mfrow=c(1,2))
      acf(diff(residuals))
      pacf(diff(residuals)) 
      sqrt(mean(residuals^2))
      #sqrt(mean(res.test^2))
      plot(diff(residuals), type='l')
      
      library("forecast")
      ar  = auto.arima(residuals,d=0)
      summary(ar)
      ar.forecast <- forecast(ar, h=48, level=c(90))
      pred.ar <- ar.forecast$mean
      lq <- matrix(ar.forecast$lower,nrow=1) ###upper bound of the 90% quantile interval
      uq <- matrix(ar.forecast$upper,nrow=1)
      
      #resi_mean = colMeans(residuals)
      #resi_std = apply(residuals,2,sd)
      # lq = qnorm(0.05,resi_mean,resi_std)
      # uq = qnorm(0.95,resi_mean,resi_std)
      # lq = matrix(rep(lq, Ns), nrow=Ns, byrow=T)
      # uq = matrix(rep(uq, Ns), nrow=Ns, byrow=T)
      
      MeanMatrix = matrix(rep(Mean, Ns), nrow=Ns, byrow=T)
      StdMatrix = matrix(rep(Std, Ns), nrow=Ns, byrow=T)
      
      forecast = (yh + matrix(pred.ar,nrow=1))*StdMatrix + MeanMatrix
      
      lowerquantile = (yh + lq)*StdMatrix + MeanMatrix
      upperquantile = (yh + uq)*StdMatrix + MeanMatrix
      
      tod = 0:47
      date = Datelist[i+1]
      
      df = data.frame(list("Date" = as.Date(date),
                           "Forecast" = c(forecast), 
                           "Load" = c(Y$DataTestY),
                           "lowerquantile" = c(lowerquantile),
                           "upperquantile" = c(upperquantile),
                           "tod"= tod))
      
      df_table = rbind(df_table,df)
    }
  }
}
#write.csv(df_table,file = file)


library(dplyr)
df_table$month =  format(as.Date(df_table$Date),"%m")
summary(df_table)



df_table$WithinRange =(df_table$Load >= df_table$lowerquantile) & (df_table$Load <= df_table$upperquantile)
CurveWithinRange = matrix(df_table$WithinRange, ncol = 48,byrow = TRUE)
df_table$CurveWithinRange = rep(rowSums(CurveWithinRange) == 48,each = 48)

CurveCR = mean(rowSums(CurveWithinRange) == 48)
PointCR = mean(df_table$WithinRange)
AvL = mean(df_table$upperquantile - df_table$lowerquantile)
cat("CurveCR:",CurveCR, "PointCR:",PointCR,"AvL:", AvL)

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

df_table %>%
  group_by(tod) %>%
  summarise("MAPE" = mean(abs((Forecast-Load)/Load)),
            "CurveCR" = mean(CurveWithinRange),
            "PointCR" = mean(WithinRange),
            "AvL" = mean(upperquantile-lowerquantile))


# table_sarx_correction2 = df_table %>%
#   group_by(tod) %>%
#   summarise("MAPE" = mean(abs((Forecast-Load)/Load)),
#             # "CurveCR" = mean(CurveWithinRangåe),
#             "PointCR" = mean(WithinRange),
#             "AvL" = mean(upperquantile-lowerquantile))
 
# DataSelect = Data
# data = PrepareData(DataSelect, feature_list,date_fin,date_fin2)
# N = dim(data$DataTrainY)[1]
# Ns = dim(data$DataTestY)[1]
# cat(N,Ns,"\n")
# Y = data
# 
# df = CurvePredictInterval(N = N, D = D, sigmaE = sigmaE, Ns = Ns, Alpha_list = Alpha_list,
#                           Ncv = Ncv, Jcv = Jcv, Nu = Nu, p = p, HO = HO, Data = "Real", Y = Y,
#                           threshold = threshold,method=method, 
#                           Dhat_s = Dhat_s)
# 
# df$weekday = "all"
# df = data.frame(df)
# df_table = rbind(df_table,df)











# file = paste("Real",feature,etaDimMulti,"result.csv",sep = '_')
# write.csv(df_table,file = file)

# 
# df_table = data.frame()
# # Real Data
# Nu = 48
# dat = read.table("../data/french_data_total.csv", header = TRUE, sep=",")
# dat$Date = as.Date(dat$dateheure) 
# dateselect = dat[dat$Instant == 8,"dateheure"]
# dat = dat[dat$Date %in% as.Date(dateselect),]
# dat = dat[1:(87401-41),]
# 
# weekdays_list = unique(dat$JourSemaine)
# 
# 
# for (method in method_list){
#   for (weekday in weekdays_list){
#     print(weekday)
#     DataSelect = dat[dat$JourSemaine == weekday,]
#     
#     DataTrain <- DataSelect[DataSelect$Annee %in% 2013:2016, ]
#     DataTest <- DataSelect[DataSelect$Annee == 2017, ]
#     
#     DataTrainY = DataTrain$Load
#     DataTrainY = matrix(DataTrainY, ncol = 48, byrow = TRUE)
#     DataTrainX = DataTrain$Load.48
#     DataTrainX = matrix(DataTrainX, ncol = 48, byrow = TRUE)
#     DataTestY = DataTest$Load
#     DataTestY = matrix(DataTestY, ncol = 48, byrow = TRUE)
#     DataTestX = DataTest$Load.48
#     DataTestX = matrix(DataTestX, ncol = 48, byrow = TRUE)
#     data = list("DataTrainY" = DataTrainY, "DataTrainX" = DataTrainX, "DataTestY" = DataTestY, "DataTestX" = DataTestX )
#     
#     N = dim(data$DataTrainY)[1]
#     Ns = dim(data$DataTestY)[1]
#     cat(N,Ns,"\n")
#     Y = data
#     
#     df = CurvePredictInterval(N = N, D = D, sigmaE = sigmaE, Ns = Ns, Alpha = Alpha,
#                               Ncv = Ncv, Jcv = Jcv, Nu = Nu, p = p, HO = HO, Data = "Real", Y = Y,threshold = threshold,method=method)
#     #df$data_type = data_type
#     df$method = method
#     df$weekday = weekday
#     df = data.frame(df)
#     df_table = rbind(df_table,df)
#   }
# }
# 

# file = paste("Real_result_old_data_weekday.csv",sep = '_')
# write.csv(df_table,file = file)
