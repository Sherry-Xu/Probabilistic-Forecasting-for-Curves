rm(list = ls())
library("glmnet")
library("polynom")
library("MASS")
source("PredictInt_F.R") 
library(foreach)
library(doParallel)

# Parallel setting
corenum = 2
myCluster <- makeCluster(corenum, # number of cores to use
                         type = "PSOCK") # type of cluster
clusterEvalQ(myCluster,source("PredictInt_F.R"))  # additional R function used inside for loop
clusterEvalQ(myCluster,library("glmnet","polynom","MASS")) #this is  R libraries used inside for loop

# Function for selecting features
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

# Setting
threshold = 0.999; p = 1; Alpha_list = c(0.60,0.01,0.10); Ncv = 600; Jcv = 6
D = NULL; sigmaE = NULL; HO = NULL
Nu = 48;

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

# Remove holidays and Christmas Break
sel <- which(Data$BH==1)
sel2 <- c(sel-48, sel, sel+48)
sel3 = c(sel2, which(Data$Summer_break!=0), which(Data$Christmas_break!=0))
Inval <- pmax(c(sel2, which(Data$Summer_break!=0), which(Data$Christmas_break!=0)),1)

length(Inval)/48
SummerBreak = Data[(Data$Summer_break != 0) & (Data$BH!=1),]
SummerBreak$Month = "SummerBreak(08)"
ChristmasBreak = Data[Data$Christmas_break != 0 & (Data$BH!=1),]
ChristmasBreak$Month = "ChristmasBreak"
Data <- Data[-Inval,]
Data[(Data$Month == '08'),'Month']  = '07'
Data = rbind(Data, SummerBreak)  #ChristmasBreak
Data = Data[order(Data$Date),]

# Datetime for training and testing dataset seperation
date_fin <- as.POSIXct(strptime("2018-12-31 23:30:00", "%Y-%m-%d %H:%M:%S"), tz="UTC")
date_fin2 <- as.POSIXct(strptime("2019-12-31 23:30:00", "%Y-%m-%d %H:%M:%S"), tz="UTC")

# selected features
feature = "17T"
feature_list = c("Load.48","Load.336", "Temp") 

# months groups
months_list = sort(unique(Data$Month))
print(months_list)
Data[Data$Month %in% c("06","07","09"), "Month"] = "06-09" 
Data[Data$Month %in% c("12","01","02"), "Month"] = "12-02"
Data[Data$Month %in% c("04","05"), "Month"] = "04-05"
months_list = sort(unique(Data$Month))
print(months_list)

# dayoftheweek groups
splitweekday = "DayofWeek"
weekdays_list = unique(Data[,splitweekday])
print(weekdays_list)

method_list = c('AIC') #,'Lasso','Ridge'
Dhat_s = "covariance" #Dhat_list[1]
Kchi2 = 5000#20000
ifplot = FALSE

# method = "AIC"
# month = "06-09"
# weekday = 'Tue-Thu'
# i = 1

df_table_all = data.frame()
file = paste("Real",Kchi2, lag, feature, Alpha_list, Dhat_s, "month_result.csv",sep = '_')

registerDoParallel(myCluster)

for (method in method_list){
  for(month in months_list){
    for (weekday in weekdays_list){
      
      tryCatch({
      print(weekday)
      print(month)
      DataSelect = Data[(Data$DayofWeek == weekday) & (Data$Month == month),] #& 
      Datelist <- sort(unique(DataSelect[(DataSelect$Date>date_fin) & (DataSelect$Date<=date_fin2),]$DateD))
      Datelist = c(as.Date(date_fin),Datelist)
      print(Datelist)
      
      
      df_table = foreach(i = 1:(length(Datelist)-1),.combine = 'rbind') %dopar% {
        #sample(1:(length(Datelist)-1), 3, replace=FALSE)
        
        print(i)
        
        date_fin_select <- as.POSIXct(strptime(paste(Datelist[i], " 23:30:00"), "%Y-%m-%d %H:%M:%S"), tz="UTC")
        date_fin_select2 <- as.POSIXct(strptime(paste(Datelist[i+1], " 23:30:00"), "%Y-%m-%d %H:%M:%S"), tz="UTC")
        
        data = PrepareData(DataSelect, feature_list, date_fin_select, date_fin_select2)
        
        N = dim(data$DataTrainY)[1]
        Ns = dim(data$DataTestY)[1]
        cat(N,Ns,"\n")
        
        Y = data
        DataType = "Real"
        date =  Datelist[i+1]
        
        df = CurvePredictInterval(N = N, D = D, sigmaE = sigmaE, Ns = Ns, Alpha_list = Alpha_list,
                                  Ncv = Ncv, Jcv = Jcv, Nu = Nu, p = p, HO = HO, DataType = DataType, Y = Y,
                                  threshold = threshold,method=method, 
                                  Dhat_s = Dhat_s, date = date, Kchi2 = Kchi2,ifplot=ifplot)
        
        
        
        df$date = Datelist[i+1]
        df$weekday = weekday
        df$Month = month
        df = data.frame(df)
        df
      }
      df_table_all = rbind(df_table_all,df_table)
      write.csv(df_table_all,file = file) 
      }, error=function(e){})
      
    }
  }
}

# Summarize the results
library(dplyr)

df_table_all$month =  format(as.Date(df_table_all$date),"%m")
summary(df_table_all)

df_table_all %>%
  group_by(month) %>%
  summarise("MAPE" = mean(MAPE.out.sample1),
            
            "Chi2.CR" = mean(condiCoverageRate.Chi2),
            "Chi2.PointCR" = mean(PointCoverageRate.Chi2),
            "Chi2.AvL" = mean(AvLChi2),
            
            "ECDF.CR" = mean(condiCoverageRate.Empirical),
            "ECDF.PointCR" = mean(PointCoverageRate.Empirical),
            "ECDF.AvL" = mean(AvLEmpirical),
            
            "ECDFB.CR" = mean(condiCoverageRate.CV.Empirical),
            "ECDFB.PointCR" = mean(PointCoverageRate.CV.Empirical),
            "ECDFB.AvL" = mean(AvLEmpiricalCV))

df_table_all %>%
  summarise("MAPE" = mean(MAPE.out.sample1),
            
            "Chi2.CR" = mean(condiCoverageRate.Chi2),
            "Chi2.PointCR" = mean(PointCoverageRate.Chi2),
            "Chi2.AvL" = mean(AvLChi2),
            
            "ECDF.CR" = mean(condiCoverageRate.Empirical),
            "ECDF.PointCR" = mean(PointCoverageRate.Empirical),
            "ECDF.AvL" = mean(AvLEmpirical),
            
            "ECDFB.CR" = mean(condiCoverageRate.CV.Empirical),
            "ECDFB.PointCR" = mean(PointCoverageRate.CV.Empirical),
            "ECDFB.AvL" = mean(AvLEmpiricalCV))
 