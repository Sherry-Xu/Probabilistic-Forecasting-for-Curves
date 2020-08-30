rm(list = ls())
library("MASS")
library("polynom")
library("glmnet")
source("PredictInt_F.R") 
source("parameters.R") 

method =  'AIC'
  
# Parallel setting
library(foreach)
library(doParallel)
myCluster <- makeCluster(corenum, # number of cores to use
                         type = "PSOCK") # type of cluster
clusterEvalQ(myCluster,source("PredictInt_F.R"))  # additional R function used inside for loop
clusterEvalQ(myCluster,source("parameters.R"))  # additional R function used inside for loop
clusterEvalQ(myCluster,library("glmnet","polynom","MASS")) #this is  R libraries used inside for loop
registerDoParallel(myCluster)


file = paste(runtime,method,dist,threshold*1000,p,N,D,sigmaE*100,HO,sparsity,Alpha_list,Kchi2, "Simu_result.csv",sep = '_')
cat(file)

# Experiments
df_table = foreach(repe=1:N_rep, .combine = 'rbind',
                   .multicombine = TRUE) %dopar% {
                     
                     set.seed(repe)
                     result = CurvePredictInterval(N = N, D = D, sigmaE = sigmaE, Ns = Ns, Alpha_list = Alpha_list,
                                                     Ncv = Ncv, Jcv = Jcv, Nu = Nu, p = p, HO = HO, DataType = DataType, Y = Y,
                                                     threshold = threshold,dist = dist,
                                                     method = method,sparsity = sparsity,Kchi2 = Kchi2,ifplot=ifplot)
  
                      df = data.frame(result)
                      df$seed = repe
                      write.table(df,file,sep=',',col.names=FALSE,append = TRUE)
                      df
                   }

write.csv(df_table,file)
summary(df_table)  
  
  