library("polynom")
library("MASS")
# library("Matrix")
# library("matrixcalc")

DataGeneration = function(N, Ns, D, Nu, sigmaE, HO, Alpha, p=1, dist = 'norm', sparsity = FALSE){
  
  v= 5
  rate_t = sqrt(sigmaE*sigmaE/(v/(v-2)))
  rate_e = sqrt(1/(sigmaE*sigmaE))
  
  BurnIn = 200
  
  U = seq(-1, 1, 2/(Nu -1)) # a grid on which values of curves are recorded
  Y = matrix(nrow = N + Ns, ncol = Nu) # Observed N curves
  XI = matrix(nrow = BurnIn + N + Ns + p, ncol = D) # random loadings of curves
  PhiT = matrix(nrow = D, ncol = Nu) # D orthonormal curves
  
  # Basis functions
  PhiT[1,] = 1/sqrt(2)
  for(i in 2:D) PhiT[i,] = cos((i-1)*3.1416*U) 
  # 1/sqrt(2), sin(k{\pi}x), cos(k{\pi}x) is a complete orthonormal basis on [-1, 1]
  
  if (HO == TRUE){
    
    if (sparsity == FALSE){
      # random loadings are AR(p*D)
      r = matrix(nrow = p*D, ncol = D)
      for (i in 1:(p*D)){
        for (j in 1:D){
          r[i,j] = runif(1, 2, 5) * (-1)**sample(c(1,2),1)
        }
      }
      parameter = matrix(0,nrow= p*D, ncol = D)
      for (i in 1:D){
        a = poly.from.roots(r[,i])
        a = (-a/a[1])[2:(p*D+1)]
        parameter[,i] = a
      }
    }else if(sparsity == 'lag'){
      #Sparsity
      r = matrix(nrow = D, ncol = D)
      for (i in 1:(D)){
        for (j in 1:D){
          r[i,j] = runif(1, 2, 5) * (-1)**sample(c(1,2),1)
        }
      }
      parameter = matrix(0,nrow= p*D, ncol = D)
      for (i in 1:D){
        a = poly.from.roots(r[,i])
        a = (-a/a[1])[2:(D+1)]
        parameter[((p-1)*D+1):(p*D),i] = a
      }
    }else if(sparsity == "diag"){
      r = matrix(nrow = p, ncol = D)
      for (i in 1:(p)){
        for (j in 1:D){
          r[i,j] = runif(1, 1.1, 5) * (-1)**sample(c(1,2),1)
        }
      }
      parameter = matrix(0,nrow= p*D, ncol = D)
      for (i in 1:D){
        a = poly.from.roots(r[,i])
        a = (-a/a[1])[2:(p+1)]
        for (j in 1:p){
        parameter[i+(j-1)*D,i] = a[j]
        }
      }
    }
    
    for (i in 1:p){
      for (j in 1:D){
        XI[i,j] = runif(1, -0.1, 0.1) * (-1)**sample(c(1,2),1)
      }
    }
    for (i in (p+1):(BurnIn+N+Ns+p)){
        # normal
      if (dist == "norm"){
        XI[i,] = matrix(t(XI[(i-p):(i-1),]),nrow=1) %*% parameter + rnorm(D,0,sigmaE)
      }else if(dist == "t"){
        # t distribution
        XI[i,] = matrix(t(XI[(i-p):(i-1),]),nrow=1) %*% parameter + rt(D, df=v)*rate_t        
      }else if(dist == "exp"){
        # exponential distribution
        XI[i,] = matrix(t(XI[(i-p):(i-1),]),nrow=1) %*% parameter + rexp(D,rate_e) - sigmaE
      }
    }
    
  }else{
    # random loadings are AR(1)
    parameter = matrix(0,nrow= D, ncol = D)
    for(i in 1:(D)) {
      a=runif(1, 0.2, 0.95)
      XI[,i]= arima.sim(n= BurnIn+N+Ns+p, list(order=c(1,0,0), ar=a*((-1)**i)), sd=sigmaE)
      parameter[i,i] = a*((-1)**i)
    }
    # # XI[,(D-1)]=rnorm(N+Ns)*sigmaE
    # # XI[,D]=rnorm(N+Ns)*sigmaE
  }
  
  # Remove BurnIn period
  XI = XI[(BurnIn+p+1):(BurnIn+N+Ns+p),]
  
  # Generate curves
  Y = XI%*%PhiT
  
  # ConfidenceInterval FromDGP
  B = 400
  
  XIE = matrix(nrow=Ns, ncol=D)
  for (i in 1:Ns){
    XIE[i,] = matrix(t(XI[(N+i-p):(N+i-1),]),nrow=1) %*% parameter
  }
  
  # normal
  if (dist == "norm"){
    residual = matrix(rnorm(B*D,0,sigmaE), nrow = B, ncol = D)
  }else if(dist == "t"){
    # t-distribution
    residual = matrix(rt(B*D, df=v)*rate_t, nrow = B, ncol = D)
  }else if(dist == "exp"){
    # exponential distribution
    residual = matrix((rexp(B*D,rate_e) - sigmaE), nrow = B, ncol = D)
  }
  
  s2Rx=vector(length=B)
  for(i in 1:B) s2Rx[i]=sum((residual[i,]/sigmaE)**2)
  
  if (dist == "norm"){
    Hu = qchisq(1-Alpha, D)
  }else{
    Nup=as.integer(B*(1-Alpha))
    Hu = sort(s2Rx, partial=Nup)[Nup]
  }
  
  # True region
  residualCI = residual[s2Rx <= Hu,]
  BCI = dim(residualCI)[1]
  
  Prange = array(dim = c(2, Ns, Nu))
  for (i in 1:Ns){
    tmp = matrix(rep(XIE[i,], BCI), nrow=BCI, byrow=T)
    XIEE = tmp + residualCI
    YE = XIEE %*%PhiT
    
    for (j in 1:Nu){
      Prange[,i,j] =  range(YE[,j])
    }
    
  }
  
  return(list("Y" = Y, "Prange" = Prange,"Parameter" = parameter))
}


PredictiveRegin = function(Dhat, Sigma, Ns, y, xih, Phi, Hu, Huempiric){
  
  xiT=y%*%Phi
  epsilon=xiT-xih # Ns x Dhat
  for(j in 1:Dhat) epsilon[,j]=epsilon[,j]/Sigma[j]
  Ii=0; Iie=0
  for(n in 1:Ns) { tmp=sum(epsilon[n,]**2)
  if(tmp>Hu) Ii=Ii+1
  if(tmp>Huempiric) Iie=Iie+1
  } 
  CoverRate=1-Ii/Ns; CoverRateE=1-Iie/Ns
  cat("Post-sample coverage rate of predictive region (Chi^2)", CoverRate, "\n")
  cat("Post-sample coverage rate of predictive region (Empricial)", CoverRateE, "\n")
  return(list("CoverageRate.Chi2"=CoverRate, "CoverageRate.Empirical"=CoverRateE))
}


ConditionalPredictiveRegionChi2 = function(Kchi2, Dhat, Sigma, Nu, Nt, Ns, y, xih, Phi, Hu, InvCov, Lmatrix, Std, ifplot){
  
  K = Kchi2
  
  Rx3 = matrix(nrow = K, ncol = Dhat)
  
  # for(j in 1:Dhat) Rx3[,j] = rnorm(K,0,Sigma[j])
  # s2Rx=vector(length= K)
  # #for(i in 1:(K)) s2Rx[i]=sum((Rx3[i,]/Sigma)**2)
  # for(i in 1:(K)) s2Rx[i]= matrix(Rx3[i,],nrow=1) %*% InvCov %*% matrix(Rx3[i,],ncol =1)
  
  for(j in 1:Dhat) Rx3[,j] = rnorm(K,0,1)
  s2Rx=vector(length= K)
  for(i in 1:(K)) s2Rx[i]= matrix(Rx3[i,],nrow=1)%*% matrix(Rx3[i,],ncol =1)
  Rx3 = t(Lmatrix %*% t(Rx3))
  
  # normRx = t(InvLmatrix %*% t(Rx3))
  # s2RxOri2=vector(length = K)
  # for(n in 1:K) s2RxOri2[n]= matrix(normRx[n,],nrow=1) %*% matrix(normRx[n,],ncol =1)
  # mean(abs(s2RxOri2 - s2Rx))
  
  Ind=rep(0,K)
  for(i in 1:K) if(s2Rx[i]>Hu) Ind[i]=1
  
  RxSub=Rx3[Ind==0,]
  RxSub = as.matrix(RxSub)
  Nr=nrow(RxSub); Prange=array(dim = c(2,Ns,Nu))
  Ii=0; Iie = 0;
  if (ifplot == TRUE){SregionChi2 = array(dim = c(Ns,Nr,Nu))}
  for(n in 1:Ns) { 
    tmp=matrix(rep(xih[n,], Nr), nrow=Nr, byrow=T)
    tmp=tmp+RxSub   # Nr x Dhat
    tmp=tmp%*%t(Phi)
    if (ifplot == TRUE){SregionChi2[n,,] = tmp}
    for(j in 1:Nu) {
      Prange[,n,j]=range(tmp[,j])
    }
    if((max(y[n,]-Prange[2,n,])>0)|(max(Prange[1,n,]-y[n,])>0)) Ii=Ii+1
    for (j in 1:Nu){
      if ((y[n,j]> Prange[2,n,j])|(Prange[1,n,j]>y[n,j])) Iie = Iie +1
    }
  }
  
  condiCoverRate = 1 - Ii/Ns
  PointwiseCoverRate = 1 - Iie/(Ns*Nu)
  
  if (Ns == 1){
    AvLChi2 = mean((Prange[2,,] - Prange[1,,]) * Std)
  }else{
    AvLChi2 = mean(colMeans(Prange[2,,] - Prange[1,,])*Std)
  }
  
  PrangeChi = Prange
  
  cat("Post-sample coverage rate of conditional predictive regions (Chi^2)", condiCoverRate, "\n")
  cat("Post-sample pointwise coverage rate of conditional predictive regions (Chi^2)", PointwiseCoverRate, "\n")
  
  
  return_list = list("CondiCoverageRate.Chi2"= condiCoverRate,
                     "PointCoverageRate.Chi2"= PointwiseCoverRate,
                     "AvLChi2" = AvLChi2,
                     "PrangeChi" = PrangeChi)
  
  if (ifplot == TRUE){
    return_list$SregionChi2 = SregionChi2
  }
  
  return(return_list)
}

BootstrapResidualCVChi2 = function(Dhat, Nu, Nt, Ns, YY, Xih, Phi, Hu, Lmatrix, Alpha){
  
  Ncv = 500
  Jcv = 5
  
  CI = matrix(rep(0, Nt*Jcv), nrow=Nt) # , ncol=2*Jcv); 
  
  K = 500 + Ncv * Jcv
  
  Rx3 = matrix(nrow = K, ncol = Dhat)
  for(j in 1:Dhat) Rx3[,j] = rnorm(K,0,1)
  s2Rx=vector(length = K)
  for(i in 1:(K)) s2Rx[i]= matrix(Rx3[i,],nrow=1)%*% matrix(Rx3[i,],ncol =1)
  Rx3 = t(Lmatrix %*% t(Rx3))
  
  for(j in 1:Jcv) {
    
    Ind=rep(0,j*Ncv); 
    for(i in 1:(j*Ncv)) { 
      if(s2Rx[i]>Hu) Ind[i]=1
    }
    
    RxS0 = Rx3[1:(j*Ncv),]
    RxS0 = as.matrix(RxS0)
    
    RxSub = RxS0[Ind==0,]; 
    RxSub = as.matrix(RxSub);
    Nr = nrow(RxSub); 
    
    Prange=array(dim = c(2,Nt,Nu))
    
    for(n in 1:Nt) { 
      tmp=matrix(rep(Xih[n,], Nr), nrow=Nr, byrow=T)
      tmp=tmp+RxSub   # Nr x Dhat
      tmp=tmp%*%t(Phi)
      for(k in 1:Nu) {
        Prange[,n,k]=range(tmp[,k])
      }
      if((max(YY[n,]-Prange[2,n,])>0)|(max(Prange[1,n,]-YY[n,])>0)) CI[n,j]=1
    }
  }
  
  CIs = vector(length = Jcv)
  colMeans(CI)
  for(j in 1:Jcv) CIs[j]=abs(sum(CI[,j])/Nt-Alpha)
  Cv1 = which.min(CIs[1:Jcv])
  print(CIs)
  Kchi = Cv1*Ncv + 500
  cat("CV: No. of  points for Chi2 region:", Kchi, "\n")
  return(list("Kchi"= Kchi))
}

ConditionalPredictiveRegionEmpirical = function(Dhat, Rx, s2RxOri, Nu, Nt, Ns, y, xih, Phi, Huempiric,InvCov, Std){
  
  # Empirical
  Ind=rep(0,Nt)
  for(n in 1:Nt) if(s2RxOri[n] > Huempiric) Ind[n] = 1
  RxSub = Rx[Ind==0,] # Nr x Dhat  residuals within 1-Alpha range
  RxSub = as.matrix(RxSub)
  Nr=nrow(RxSub); Prange=array(dim=c(2,Ns,Nu))
  Ii=0; Iie = 0
  SregionEmpirical = array(dim=c(Ns,Nr,Nu))
  for(n in 1:Ns) { 
    tmp = matrix(rep(xih[n,], Nr), nrow=Nr, byrow=T)
    tmp = tmp + RxSub   # Nr x Dhat
    tmp = tmp%*%t(Phi)
    SregionEmpirical[n,,] = tmp
    for(j in 1:Nu) {
      Prange[,n,j]=range(tmp[,j])
    }
    if((max(y[n,]-Prange[2,n,])>0)|(max(Prange[1,n,]-y[n,])>0)) Ii=Ii+1
    for (j in 1:Nu){
      if ((y[n,j] > Prange[2,n,j])|(Prange[1,n,j]>y[n,j])) Iie = Iie +1
    }
  }
  
  condiCoverRateE=1-Ii/Ns
  PointwiseCoverRate = 1 - Iie/(Ns*Nu)
  
  if (Ns == 1){
    AvLEmpirical = mean((Prange[2,,] - Prange[1,,])*Std)
  }else{
    AvLEmpirical = mean(colMeans(Prange[2,,] - Prange[1,,])*Std)
  }
  
  PrangeEmpirical = Prange
  
  cat("Post-sample coverage rate of conditional predictive regions (Empirical)", condiCoverRateE, "\n")
  cat("Post-sample pointwise coverage rate of conditional predictive regions (Empirical)", PointwiseCoverRate, "\n")
  
  return(list("CondiCoverageRate.Empirical" = condiCoverRateE,
              "PointCoverageRate.Empirical" = PointwiseCoverRate,
              "AvLEmpirical" = AvLEmpirical, 
              "PrangeEmpirical" = PrangeEmpirical,
              "SregionEmpirical"= SregionEmpirical))
}

BootstrapResidualCVEmpirical = function(Dhat,Sigma, Rx, Jcv, Ncv, Nt, Nu, YY, Xih, Phi, Huempiric, InvCov, Alpha){
  
  CI=matrix(rep(0, Nt*Jcv), nrow=Nt) # , ncol=2*Jcv); 
  
  for(n in 1:Nt) {   
    #Bootstrap K, Points
    Rx1=Rx[-n,]; K=(Jcv-1)*Ncv
    Rx1 = as.matrix(Rx1)
    Rx2 = matrix(nrow = K, ncol = Dhat)
    for(j in 1:Dhat) Rx2[,j] = sample(Rx1[,j], K, replace=T) # Wrong?
    tt=data.frame(t(Rx1), t(Rx2))
    Rx3=as.matrix(t(tt))
    
    s2Rx = vector(length= dim(Rx3)[1])
    #for(i in 1:dim(Rx3)[1]) s2Rx[i]=sum((Rx3[i,]/Sigma)**2)
    for(i in 1:dim(Rx3)[1]) s2Rx[i]= matrix(Rx3[i,],nrow=1) %*% InvCov %*% matrix(Rx3[i,],ncol =1)
    
    for(j in 1:Jcv) {
      # Nup=as.integer((Nt-1+(j-1)*Ncv)*(1-Alpha))
      # Huempiric=sort(s2Rx[1:(Nt-1+(j-1)*Ncv)], partial=Nup)[Nup]
      
      Ind=rep(0,(Nt-1+(j-1)*Ncv)); 
      for(i in 1:(Nt-1+(j-1)*Ncv)) { 
        if(s2Rx[i]>Huempiric) Ind[i]=1
      }
      
      RxS0 = Rx3[1:(Nt-1+(j-1)*Ncv),]
      RxS0 = as.matrix(RxS0)
      
      RxSub = RxS0[Ind==0,]; 
      RxSub = as.matrix(RxSub);
      Nr = nrow(RxSub); 
      Prange = matrix(nrow=Nu, ncol=2); 
      tmp = matrix(rep(Xih[n,], Nr), nrow=Nr, byrow=T)
      tmp = tmp+RxSub   # Nr x Dhat
      tmp = tmp%*%t(Phi);
      for(i in 1:Nu) {Prange[i,] = range(tmp[,i]) }
      if((max(YY[n,]-Prange[,2])>0)|(max(Prange[,1]-YY[n,])>0)) CI[n,j]=1
    }
  }
  
  CIs = vector(length = Jcv)
  for(j in 1:Jcv) CIs[j]=abs(sum(CI[,j])/Nt-Alpha)
  Cv1=which.min(CIs[1:Jcv])-1
  print(CIs)
  cat("CV: No. of multiple points to add for Empirical region:", Cv1, "\n")
  return(list("Cv1"=Cv1))
}

ConditionalPredictiveReginCVEmpirical = function(Dhat, Cv1, Sigma, Rx, Ncv, Nu, Nt, Ns, y, xih, Phi, Huempiric, InvCov, Alpha,Std){
  
  K=Cv1*Ncv
  
  # Bootstrap K, points
  if(Cv1>0){
    Rx2 = matrix(nrow=K, ncol=Dhat)
    for(j in 1:Dhat) Rx2[,j] = sample(Rx[,j], K, replace=T)
    tt = data.frame(t(Rx), t(Rx2))
    Rx3 = as.matrix(t(tt))
  } else Rx3 = Rx
  
  s2Rx=vector(length=Nt+K)
  #for(i in 1:(Nt+K)) s2Rx[i]=sum((Rx3[i,]/Sigma)**2)
  for(i in 1:(Nt+K)) s2Rx[i]= matrix(Rx3[i,],nrow=1) %*% InvCov %*% matrix(Rx3[i,],ncol =1)
  # Nup=as.integer((Nt+K)*(1-Alpha))
  # Huempiric=sort(s2Rx, partial=Nup)[Nup]
  
  Ind=rep(0,(Nt+K))
  for(i in 1:(Nt+K)) if(s2Rx[i]>Huempiric) Ind[i]=1
  RxSub=Rx3[Ind==0,]
  RxSub = as.matrix(RxSub)
  Nr=nrow(RxSub); Prange= array(dim = c(2,Ns,Nu))
  Ii=0;Iie=0
  SregionEmpiricalCV = array(dim = c(Ns,Nr,Nu))
  for(n in 1:Ns) { 
    tmp=matrix(rep(xih[n,], Nr), nrow=Nr, byrow=T)
    tmp=tmp+RxSub   # Nr x Dhat
    tmp=tmp%*%t(Phi)
    SregionEmpiricalCV[n,,] = tmp
    for(j in 1:Nu){
      Prange[,n,j]=range(tmp[,j])
    } 
    if((max(y[n,]-Prange[2,n,])>0)|(max(Prange[1,n,]-y[n,])>0)) Ii=Ii+1
    for (j in 1:Nu){
      if ((y[n,j] > Prange[2,n,j])|(Prange[1,n,j] > y[n,j])) Iie = Iie +1
    }
  }
  
  condiCoverRateCV=1-Ii/Ns
  PointwiseCoverRateCV = 1 - Iie/(Ns*Nu)
  
  if (Ns == 1){
    AvLEmpiricalCV = mean((Prange[2,,] - Prange[1,,])*Std)
  }else{
    AvLEmpiricalCV = mean(colMeans(Prange[2,,] - Prange[1,,])*Std)
  }
  
  PrangeEmpirical = Prange
  
  cat("Post-sample coverage rate of bootstrapped conditional predictive regions (Empirical)", condiCoverRateCV, "\n\n\n")
  cat("Post-sample pointwise coverage rate of bootstrapped conditional predictive regions (Empirical)", PointwiseCoverRateCV, "\n\n\n")
  
  return(list("CondiCoverageRateCV.Empirical"= condiCoverRateCV,
              "PointCoverageRateCV.Empirical"= PointwiseCoverRateCV,
              "AvLEmpiricalCV" = AvLEmpiricalCV,
              "PrangeEmpiricalCV" = PrangeEmpirical,
              "Huempiric" = Huempiric,
              "SregionEmpiricalCV" = SregionEmpiricalCV))
}


CurvePredictInterval = function(N, D, sigmaE, Ns, Alpha_list = 0.1, Ncv = 200, Jcv=6,
                                Nu = 51, p = 1, HO = FALSE, DataType = "Sim", Y = NULL,threshold = 0.95,
                                dist = "norm", method = "AIC", Dhat_s = NULL,
                                date = NULL, sparsity = FALSE, ifplot = FALSE, Kchi2 = 10000) {
  # N: No of observed curves
  # Ns: test sample size
  # Nu: No of the point over the interval [-1, 1]
  # Alpha: predictive region with converage probability 1-Alpha 
  # Ncv, Jcv: Choosing the number of bootstrap samples to be N + j x Ncv for j=0,1, ..., (Jcv-1)
  # p: FAR(p) model
  # Simulated curves are defined on the interval [-1, 1]
  # D: curve dimension
  # sigmaE: STD of noise in linear regression
  # HO: Higher order of xi_tj on eta_tl
  # Y: Real Data
  
  start.time <- Sys.time()
  
  # Part I: Generate Curve Series
  if (DataType == "Sim"){
    
    DGP = DataGeneration(N, Ns, D, Nu, sigmaE,HO, 0.1, p, dist,sparsity)
    
    Y = DGP$Y
    TruePrange = DGP$Prange
    
    Mean = colMeans(Y[1:N,])

    TruePrange[1,,] = TruePrange[1,,] - matrix(rep(Mean, Ns), nrow=Ns, byrow=T) 
    TruePrange[2,,] = TruePrange[2,,] - matrix(rep(Mean, Ns), nrow=Ns, byrow=T) 
    TrueAvL = mean(TruePrange[2,,] - TruePrange[1,,])
    
    Ydm = Y - matrix(rep(Mean, N+Ns), nrow=N+Ns, byrow=T) # de-meaned Y
    
    # training data
    YY = Ydm[(p+1):N,]
    XX = matrix(nrow = N-p, ncol = p*Nu)
    for (i in 1:p){
      XX[,((i-1)*Nu+1):(i*Nu)] = Ydm[i:(N-p+i-1),]
    }
    
    # testing data
    y= Ydm[(N+1):(N+Ns),]
    x = matrix(nrow = Ns, ncol = p*Nu)
    for (i in 1:p){
      x[,((i-1)*Nu+1):(i*Nu)] = Ydm[(N-p+i):(N+Ns-p+i-1),]
    } 
    
    Std = rep(1,Nu)
    Stdx = rep(1, dim(XX)[2])
    
  }else{
    
    YY = Y$DataTrainY
    XX = Y$DataTrainX
    y = Y$DataTestY
    x = Y$DataTestX
    
    Mean = colMeans(YY)
    Meanx = colMeans(XX)
    
    YY = (YY - matrix(rep(Mean, N), nrow=N, byrow=T))
    XX = (XX - matrix(rep(Meanx, N), nrow=N, byrow=T))
    y = (y - matrix(rep(Mean, Ns), nrow=Ns, byrow=T))
    x = (x - matrix(rep(Meanx, Ns), nrow=Ns, byrow=T))
    
    Std = apply(YY, 2, sd)
    Stdx = apply(XX, 2, sd)
    
    YY = YY /matrix(rep(Std, N), nrow=N, byrow=T)
    XX = XX/matrix(rep(Stdx, N), nrow=N, byrow=T)
    y = y /matrix(rep(Std, Ns), nrow=Ns, byrow=T)
    x = x/matrix(rep(Stdx, Ns), nrow=Ns, byrow=T)
  }
  
  Nt = dim(YY)[1] #N-p
  
  # Part II: SVD 
  # SVD for Cov(YY, XX)
  Syx =cov(YY,XX)  # sample Cov(Y, X)
  #Syx = cor(YY,XX)
  TT=svd(Syx)
  la=TT$d
  print(la)
  d0 = min(Nu,round(Nt/2))
  la1 = (la[2:Nu]/la[1:(Nu-1)])[1:d0]
  Dhat1 = which.min(la1)
  cat("Correlation dimension choice 1 :", Dhat1, "\n");
  la2 = cumsum(la)/sum(la)
  Dhat2 = sum(la2 < threshold) + 1
  cat("Correlation dimension choice 1 :", Dhat2, "\n");
  Dhat = max(Dhat1, Dhat2) #min(Nt/2,Dhat2) #max(Dhat1, Dhat2) #Dhat2 
  #Dhat = Dhat_s
  cat("Correlation dimension:", Dhat, "\n");
  if (Dhat >= Nt){
    Dhat = Nt - 1
  }
  # if (Dhat >= Nt/2){
  #   Dhat = round(Nt/2)
  # }
  # if (Dhat == 1){
  #   Dhat = 2
  # }
  cat("Correlation dimension:", Dhat, "\n");
  
  etaDhat = dim(TT$v)[2]#etaDimMulti * Dhat
  if (DataType == "Real"){
    etaDhat = min(round(1/2*Nt),etaDhat)
    if (etaDhat < Dhat){
      etaDhat = Dhat
    }
  }else{
    etaDhat = 10
  }
  cat("Eta dimension:",etaDhat)
  
  Phi =TT$u[,1:Dhat]; Psi = TT$v[,1:etaDhat] # Two sets of estimated orthonormal functions: Nu x Dhat
  Phi = as.matrix(Phi); Psi = as.matrix(Psi)
  
  # plots of the components
  # par(mfrow=c(2,Dhat), mar=c(3,4,1,1))
  # for(k in 1:Dhat)  plot(Phi[,k], type="l", xlab="", ylab=expression(phi))
  # for(k in 1:Dhat)  plot(Psi[,k], type="l", xlab="", ylab=expression(psi))
  
  Xi = YY%*%Phi  # Estimated loading for Y: (N-1) x Dhat 
  Eta = XX%*%Psi # Estimated loading for X: (N-1) x Dhat
  
  # Part III: linear regression of \xi_tj on \eta_tl

  Beta = matrix(0, nrow = etaDhat, ncol = Dhat)
  Sigma = vector(length = Dhat) #standard deviation
  Rx = matrix(nrow = Nt, ncol = Dhat) #residuals
  squareR = vector(length=Dhat)
  
  if (method == "AIC"){
    print("AIC")
    #Beta = data.frame(Beta)
    for(j in 1:Dhat) {
      variables = names(data.frame(Eta))
      model0 = lm(Xi[,j] ~ 0 + .,data = data.frame(Eta))
      lowermodel =  paste("~",paste(c(0,variables[j]), collapse="+"))
      #print(lowermodel)
      model1 = stepAIC(model0, direction = 'both',trace = 0, scope= list(lower = lowermodel))
      summary(model1)
      coef = model1$coefficients
      for (name in variables[!is.element(variables, names(coef))]){
        coef[name] = 0
      }
      coef = coef[variables]
      Beta[,j] = coef # Estimated slope of j-th regression
      Sigma[j]=(summary(model1)$sigma)    # STD of residuals for j-th regression
      squareR[j]=summary(model1)$r.squared   # squared regression correlation
      Rx[,j]=model1$residuals    # residuals of j-th regression
    
      # Residual Checking
      # standardResi =  model1$residuals/(summary(model1)$sigma)
      # par(mfrow=c(2,2))
      # plot(standardResi)
      # acf(standardResi)
      # hist(standardResi)
      # qqnorm(standardResi); qqline(standardResi)
      # Boxtest = Box.test(standardResi,lag=5,type="Ljung")
      # Shapirotest = shapiro.test(standardResi)
      # if (Boxtest$p.value >= 0.05) {correlation = "not correlated"
      # }else {correlation = "correlated"
      # }
      # if (Shapirotest$p.value >= 0.05){normality = "normal"
      # }else {normality = "not normal"
      # }
      # mtext(side  = 1 , paste("Dimension:", j, ",", correlation,",", normality, "\n"),outer=TRUE)
      
    }
  }else if (method == "Lasso"){
    print("Lasso")
    for(j in 1:Dhat) {
      cv = cv.glmnet(Eta, Xi[,j], alpha = 1,intercept= FALSE)
      #print(cv$lambda)
      cat(cv$lambda.min,length(cv$lambda),which(cv$lambda == cv$lambda.min),"\n")
      model1 <- glmnet(Eta, Xi[,j], alpha = 1, lambda = cv$lambda.min,intercept= FALSE)
      coeffi = as.vector(coef(model1))[2:(etaDhat+1)]
      residual = Xi[,j] - predict(model1,Eta)
      #print(coeffi)
      #print(model1$df)
      Beta[,j] = coeffi # Estimated slope of j-th regression
      Sigma[j]= sqrt(sum(residual ** 2)/ (length(Xi[,j]) - model1$df))    # STD of residuals for j-th regression
      squareR[j]= model1$dev.ratio   # squared regression correlation
      Rx[,j]= residual   # residuals of j-th regression
    }
  }else if (method == "Ridge"){
    print("Ridge")
    for(j in 1:Dhat) {
      cv = cv.glmnet(Eta, Xi[,j], alpha = 0,intercept= FALSE)
      #print(cv$lambda)
      #print(cv$lambda.min)
      cat(cv$lambda.min,length(cv$lambda),which(cv$lambda == cv$lambda.min),"\n")
      model1 <- glmnet(Eta, Xi[,j], alpha = 0, lambda = cv$lambda.min,intercept= FALSE)
      coeffi = as.vector(coef(model1))[2:(etaDhat+1)]
      residual = Xi[,j] - predict(model1,Eta)
      #print(coeffi)
      #print(model1$df)
      Beta[,j] = coeffi # Estimated slope of j-th regression
      Sigma[j]=sqrt(sum(residual ** 2)/ (length(Xi[,j]) - model1$df))   # STD of residuals for j-th regression
      squareR[j]= model1$dev.ratio   # squared regression correlation
      Rx[,j]= residual   # residuals of j-th regression
    }
  }else{
    print("Choose method to be AIC, Lasso or Ridge")
  }
  
  print(dim(Beta))
  sparsityall = sum(abs(Beta) > 0.0001)
  sparsitydiag = sum(abs(diag(Beta)) > 0.0001)
  sparsityrate = sparsityall/(Dhat*etaDhat)
  sparsitydiagrate = sparsitydiag/Dhat
  cat("Sparsity:",sparsityrate,",",sparsitydiagrate)
  cat("Residual STD and Squared regression correlation:", "\n")
  for(j in 1:Dhat) cat(Sigma[j], squareR[j], "\n")
  
  # Part IV: calculate MAD for prediction
  # Calculate in-sample predictive errors
  Beta = data.matrix(Beta)
  Xih = Eta%*%Beta
  YYh = Xih%*%t(Phi)
  
  MeanMatrix = matrix(rep(Mean, Nt), nrow=Nt, byrow=T)
  StdMatrix = matrix(rep(Std, Nt), nrow=Nt, byrow=T)
  
  MADinsample = mean(abs(YY*StdMatrix-YYh*StdMatrix))
  RMSEinsample = sqrt(mean((YY*StdMatrix-YYh*StdMatrix)**2))
  MAPEinsample = abs((YYh*StdMatrix+MeanMatrix)/(YY*StdMatrix+MeanMatrix)-1)
  MAPEinsample[MAPEinsample == Inf] = NA
  MAPEinsample1 = mean(MAPEinsample,na.rm=TRUE)
  MAPEinsample2 = median(MAPEinsample,na.rm=TRUE)
  cat("In-sample MAD, RMSE, MAPE1, MAPE2:", MADinsample,RMSEinsample,MAPEinsample1,MAPEinsample2,"\n")
  
  # Calculate out-sample predictive errors
  eta = x%*%Psi
  xih = eta%*%Beta
  yh =  xih%*%t(Phi)
  
  MeanMatrix = matrix(rep(Mean, Ns), nrow=Ns, byrow=T)
  StdMatrix = matrix(rep(Std, Ns), nrow=Ns, byrow=T)
  
  MAD = mean(abs(y*StdMatrix-yh*StdMatrix))
  RMSE = sqrt(mean((y*StdMatrix-yh*StdMatrix)**2))
  MAPE = abs((yh*StdMatrix+MeanMatrix)/(y*StdMatrix+MeanMatrix)-1)
  MAPE[MAPE == Inf] = NA
  MAPE1 = mean(MAPE,na.rm=TRUE)
  MAPE2 = median(MAPE,na.rm=TRUE)
  cat("Post MAD, RMSE, MAPE1, MAPE2:", MAD,RMSE,MAPE1,MAPE2,"\n")
  
  Alpha_len = length(Alpha_list)
  par(mfrow=c(Alpha_len,3))
  
  if (ifplot == TRUE){
    PrangeChi2_all = array(dim=c(2,Ns,Nu,Alpha_len))
    PrangeEmpirical_all = array(dim=c(2,Ns,Nu,Alpha_len))
    PrangeEmpiricalCV_all = array(dim=c(2,Ns,Nu,Alpha_len))
    SregionChi2 = list()
    SregionEmpiricalCV = list()
  }
  
  for (l in 1:Alpha_len){
    
    Alpha = Alpha_list[l]
    
    ## Part V: calculate predictive region
    Hu = qchisq(1-Alpha, Dhat)
    print(Hu)
    
    dim(Beta)
    
    # Covariance matrix df = (i+j)/2  (1)
    # df_free = vector(length = Dhat)
    # for (j in 1:Dhat){df_free[j] = sum(Beta[,j] != 0)}
    # covariance = matrix(0,nrow=Dhat,ncol=Dhat)
    # for (i in 1:Dhat){
    #   for (j in 1:Dhat){
    #     covariance[i,j] = sum(Rx[,i] * Rx[,j])/(Nt - (df_free[i]+df_free[j])/2)
    #   }
    # }
    
    # Covariance matrix df = i*j  (2)
    # df_free = matrix(0,nrow=Dhat,ncol=Dhat)
    # covariance = matrix(0,nrow=Dhat,ncol=Dhat)
    # for (i in 1:Dhat){
    #   for (j in 1:Dhat){
    #     df_free[i,j] = sum((Beta[,i] != 0)&(Beta[,j] != 0))
    #     covariance[i,j] = sum(Rx[,i] * Rx[,j])/(Nt - df_free[i,j])
    #   }
    # }

    #Covariance matrix df = (i+j)  (3)
    df_free = matrix(0,nrow=Dhat,ncol=Dhat)
    covariance = matrix(0,nrow=Dhat,ncol=Dhat)
    for (i in 1:Dhat){
      for (j in 1:Dhat){
        df_free[i,j] = sum((Beta[,i] != 0)|(Beta[,j] != 0))
        covariance[i,j] = sum(Rx[,i] * Rx[,j])/(Nt - df_free[i,j])
      }
    }
    
    #Covariance matrix  (4)
    # df_free = etaDhat
    # for (i in 1:etaDhat){
    #   if (sum(Beta[i,] == 0) == Dhat) df_free = df_free -1
    # }
    # covariance = matrix(0,nrow=Dhat,ncol=Dhat)
    # for (i in 1:Dhat){
    #   for (j in 1:Dhat){
    #     covariance[i,j] = sum(Rx[,i] * Rx[,j])/(Nt - df_free)
    #   }
    # }
    
    # Mixed covariance matrix
    #covariance = cov(Rx,Rx)
    # diag(covariance) = Sigma*Sigma
    
    # Diagnal variance Matrix (0)
    # covariance = matrix(0,nrow=Dhat,ncol=Dhat)
    # diag(covariance) = Sigma*Sigma

    InvCov = solve(covariance)
    s2RxOri=vector(length = Nt)
    for(n in 1:Nt) s2RxOri[n]= matrix(Rx[n,],nrow=1) %*% InvCov %*% matrix(Rx[n,],ncol =1)
    Nup=as.integer(Nt*(1-Alpha));
    Huempiric=sort(s2RxOri, partial=Nup)[Nup]
    print(Huempiric)
    
    print(eigen(covariance)$values) 
    if (min(eigen(covariance)$values) < 0){
      print("make covariance PD")
      m_threshold = 0.000001
      m_add = (m_threshold - min(eigen(covariance)$values))
      covariance = covariance + m_add * diag(Dhat)
      print(eigen(covariance)$values)
    }else{
      m_add = 0
    }
    
    Lmatrix = t(chol(covariance))
    InvLmatrix = solve(Lmatrix)
   
    # normRx = t(InvLmatrix %*% t(Rx))
    # s2RxOri2=vector(length = Nt)
    # for(n in 1:Nt) s2RxOri2[n]= matrix(normRx[n,],nrow=1) %*% matrix(normRx[n,],ncol =1)
    
    # s2RxOri=vector(length = Nt)
    # for(n in 1:Nt) s2RxOri[n]=sum((Rx[n,]/Sigma)**2)
    # Nup=as.integer(Nt*(1-Alpha));
    # Huempiric=sort(s2RxOri, partial=Nup)[Nup]
    # print(Huempiric)
    
    # CR = PredictiveRegin(Dhat, Sigma, Ns, y, xih, Phi, Hu, Huempiric)
    # CoverRate = CR$CoverageRate.Chi2
    # CoverRateE = CR$CoverageRate.Empirical
    
    # if (DataType == "Sim"){
    #   Kchi2 = BootstrapResidualCVChi2(Dhat, Nu, Nt, Ns, YY, Xih, Phi, Hu, Lmatrix, Alpha)
    #   Kchi2 = Kchi2$Kchi
    # }
    
    # Part VI: calculate conditional predictive region
    CPRChi2 = ConditionalPredictiveRegionChi2(Kchi2,Dhat, Sigma, Nu, Nt, Ns, y, xih, Phi, Hu, InvCov, Lmatrix, Std, ifplot)
    condiCoverRate = CPRChi2$CondiCoverageRate.Chi2
    PointCoverRate = CPRChi2$PointCoverageRate.Chi2
    AvLChi2 = CPRChi2$AvLChi2
    PrangeChi = CPRChi2$PrangeChi
    
    CPREmpirical =  ConditionalPredictiveRegionEmpirical(Dhat, Rx, s2RxOri, Nu, Nt, Ns, y, xih, Phi, Huempiric, InvCov,Std)
    condiCoverRateE = CPREmpirical$CondiCoverageRate.Empirical
    PointCoverRateE = CPREmpirical$PointCoverageRate.Empirical
    AvLEmpirical = CPREmpirical$AvLEmpirical
    PrangeEmpirical = CPREmpirical$PrangeEmpirical
    
    # Part VII: CV choose number of points in (2.17) 
    if (DataType == "Sim"){
      CvK = BootstrapResidualCVEmpirical(Dhat,Sigma, Rx, Jcv, Ncv, Nt, Nu, YY, Xih, Phi, Huempiric, InvCov, Alpha)
      Cv1 = CvK$Cv1
    }else{
    Cv1 = 5
    }
    
    # Part VIII: conditional predictive region determined by bootstrap (CV)
    CPRCV = ConditionalPredictiveReginCVEmpirical(Dhat, Cv1, Sigma, Rx, Ncv, Nu, Nt, Ns, y, xih, Phi, Huempiric, InvCov, Alpha,Std)
    condiCoverRateCVE = CPRCV$CondiCoverageRateCV.Empirical
    PointCoverRateCVE = CPRCV$PointCoverageRateCV.Empirical
    AvLEmpiricalCV = CPRCV$AvLEmpiricalCV
    PrangeEmpiricalCV = CPRCV$PrangeEmpiricalCV
    HuempiricCV = CPRCV$Huempiric

    if (DataType == "Sim"){
      
      a0 = pmin(PrangeChi[2,,], TruePrange[2,,])
      b0 = pmax(PrangeChi[1,,], TruePrange[1,,])
      c0 = a0-b0
      c0[c0 < 0] = 0
      OverlapChi = mean(c0)
      
      a0 = pmin(PrangeEmpirical[2,,], TruePrange[2,,])
      b0 = pmax(PrangeEmpirical[1,,], TruePrange[1,,])
      c0 = a0-b0
      c0[c0 < 0] = 0
      OverlapEmpirical = mean(c0)
      
      a0 = pmin(PrangeEmpiricalCV[2,,], TruePrange[2,,])
      b0 = pmax(PrangeEmpiricalCV[1,,], TruePrange[1,,])
      c0 = a0-b0
      c0[c0 < 0] = 0
      OverlapEmpiricalCV = mean(c0)
    }
    
    if (ifplot == TRUE){
        SregionChi2[[l]] =  CPRChi2$SregionChi2
        SregionEmpiricalCV[[l]] = CPRCV$SregionEmpiricalCV
        PrangeChi2_all[,,,l] = PrangeChi
        PrangeEmpirical_all[,,,l] = PrangeEmpirical
        PrangeEmpiricalCV_all[,,,l] = PrangeEmpiricalCV
    }
    
  }
  
  cat(condiCoverRate,condiCoverRateE,condiCoverRateCVE,'\n')
  cat(PointCoverRate,PointCoverRateE,PointCoverRateCVE,'\n')
  cat(AvLChi2,AvLEmpirical,AvLEmpiricalCV,'\n')
  
  if (ifplot == TRUE){
    if (DataType == "Real"){
      # Out of sample forecast plot
      n_list = c(1:Ns)
      print(n_list)
  
      
      # filename = paste(date,weekdays(date),"Chi2.png",sep='_')
      # filename = paste("images/",filename)
      # png(filename= filename,width = 5, height = 5, units = 'in',res=300)
      filename = paste(date,weekdays(date),"Chi2.pdf",sep='_')
      filename = paste("images/",filename)
      pdf(file= filename,width = 5, height = 5)
      par(mfrow=c(1,1),mai=c(0.8,0.8,0.6,0.2))
      for (n in n_list){
        ConfidenceIntervalPlot3(SregionChi2, PrangeChi2_all, y, n, Nu, Alpha_list, Mean = Mean, Std = Std, paste(date,weekdays(date)), yhat = yh, PointCR = PointCoverRate)
      }
      dev.off()
  
      # # filename = paste(date,weekdays(date),"ECDF(Bootstrap).png",sep='_')
      # # filename = paste("images/",filename)
      # # png(filename= filename,width = 5, height = 5, units = 'in',res=300)
      # filename = paste(date,weekdays(date),"ECDF(Bootstrap).pdf",sep='_')
      # filename = paste("images/",filename)
      # pdf(file= filename,width = 5, height = 5)     
      # for (n in n_list){
      #   ConfidenceIntervalPlot3(SregionEmpiricalCV,PrangeEmpiricalCV_all,y, n, Nu, Alpha_list, Mean, Std = Std, paste(date,weekdays(date),"ECDF-B"), yhat = yh, PointCR = PointCoverRateCVE)
      # }
      # dev.off()
    }else{
      
      n_list = sample(1:Ns,15)
      n_list

      
      for (n in n_list){
        filename = paste("Simulation_Chi2",n,".png",sep='_')
        filename = paste("images/",filename)
        png(filename= filename,width = 5, height = 5, units = 'in',res=200)
        par(mfrow=c(1,1),mai=c(0.5,0.5,0.3,0.2))
        ConfidenceIntervalPlot4(TruePrange,SregionChi2, PrangeChi2_all, y, n, Nu, Alpha_list, Mean = Mean, Std = Std, "Chi-square", yhat = yh)
        dev.off()
      }
      
      
      # # filename = paste("Simulation_Chi2_ECDFB",n,".png",sep='_')
      # # filename = paste("images/",filename)
      # # png(filename= filename,width = 5, height = 5, units = 'in',res=100)
      # # par(mfrow=c(2,5))
      # for (n in n_list){
      #   filename = paste("Simulation_Chi2_ECDFB",n,".png",sep='_')
      #   filename = paste("images/",filename)
      #   png(filename= filename,width = 5, height = 5, units = 'in',res=200)
      #   par(mfrow=c(1,1),mai=c(0.5,0.5,0.3,0.2))
      #   ConfidenceIntervalPlot4(TruePrange,SregionEmpiricalCV,PrangeEmpiricalCV_all,y, n, Nu, Alpha_list, Mean, Std = Std, "ECDF-B", yhat = yh)
      #   dev.off()
      #   }
      
    }
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  # result1 = list("Corr.dim1"=Dhat1, "Corr.dim2"=Dhat2, "Corr.dim"= Dhat, "Phi"=Phi, "Psi"=Psi, "Mean"=Mean, "Resid.STD"=Sigma, "RegressR2"=squareR, "MAD.in.sample"=MADinsample, "MAD.out.sample"=MAD, "CoverageRate.Chi2"=CoverRate, "condiCoverageRate.Chi2"=condiCoverRate, "condiCoverageRate-CV.Chi2"=condiCoverRateCV, "CoverageRate.Empirical"=CoverRateE, "condiCoverageRate.Empirical"=condiCoverRateE, "condiCoverageRate-CV.Empirical"=condiCoverRateCVE, "CV-AddPoints.Chi2"=Cv1*Ncv, "CV-AddPoints.Empirical"=Cv2*Ncv)
  result = list("Method" = method, "N" = N, "Ns" = Ns,
                "Corr.dim1"= Dhat1, "Corr.dim2"= Dhat2, "Corr.dim"= Dhat, "etaDhat" = etaDhat,
                "MAD.in.sample" = MADinsample, "MAD.out.sample" = MAD, 
                "MAPE.in.sample1" = MAPEinsample1, "MAPE.out.sample1" = MAPE1, 
                "MAPE.in.sample2" = MAPEinsample2, "MAPE.out.sample2" = MAPE2, 
                "RMSE.in.sample"  = RMSEinsample, "RMSE.out.sample" = RMSE, 
                "condiCoverageRate.Chi2" = condiCoverRate,
                "PointCoverageRate.Chi2" = PointCoverRate,
                "condiCoverageRate.Empirical" = condiCoverRateE, 
                "PointCoverageRate.Empirical" = PointCoverRateE, 
                "condiCoverageRate-CV.Empirical" = condiCoverRateCVE,
                "PointCoverageRate-CV.Empirical" = PointCoverRateCVE, 
                "AvLChi2" = AvLChi2,
                "AvLEmpirical" = AvLEmpirical,
                "AvLEmpiricalCV" = AvLEmpiricalCV,
                "CV-AddPoints.Empirical"=Cv1*Ncv, 
                "TimeTaken" = time.taken, 
                "Sparsity" = sparsityall, "SparsityRate" = sparsityrate,
                "SparsityDiag" = sparsitydiag, "SparsityDiagRate" = sparsitydiagrate,
                "Hu" = Hu, "Huempiric" = Huempiric, "HuempiricCV" = HuempiricCV,
                "mean" = mean(Mean),"m_add"= m_add,"Kchi2" = Kchi2)
  
  if (DataType == "Sim"){
    result_add = list("TrueAvL" = TrueAvL,
                      "OverlapChi" = OverlapChi,
                      "OverlapEmpirical" = OverlapEmpirical,
                      "OverlapEmpiricalCV" = OverlapEmpiricalCV)
    result = c(result,result_add)
  }
  
  return(result)
}  
   

ConfidenceIntervalPlot3 =  function(S_list, Prange_all, y_all, n, Nu, Alpha_list, Mean= NULL, Std = NULL,
                                    title_name = "PredictiveRegion", yhat = NULL,PointCR = NULL){
  
  #N = dim(Prange_all)[2]
  
  Alpha_len = length(Alpha_list)
  
  Prange = array(dim=c(dim(Prange_all)[1],dim(Prange_all)[3],dim(Prange_all)[4]))
  Prange[,,] = Prange_all[,n,,]
  
  y = y_all[n,]
  yhat = yhat[n,]
  
  if (!is.null(Mean)){
    Mean2 = matrix(rep(Mean,2),nrow = 2, byrow = TRUE)
    Std2 = matrix(rep(Std,2),nrow = 2, byrow = TRUE)
    for (l in 1:Alpha_len){
      Prange[,,l] = Prange[,,l] * Std2 + Mean2
    }
    y =  y*Std + Mean
    yhat = yhat * Std + Mean
  }
  
  if((max(y-Prange[2,,2])>0)|(max(Prange[1,,2]-y)>0)){
    Status = "Out of 90% Predictive Region,"
    Status = paste(Status, "pointwiseCR:", round(PointCR,3))
  }else{
    Status = "Within 90% Predictive Region"
  }
  
  lb = min(Prange)
  ub = max(Prange)
  
  xind = (0:(Nu-1))/2
  xind2 = c(xind[25:48],xind[1:24])
  #par(mfrow= c(1,1))
  plot(xind, y, xaxt="n",xlab="hour",ylab="load", col='black', ylim=c(lb, ub), type='l',lwd=3)
  
  
  l = 3
  Alpha = Alpha_list[l]
  zmax = Prange[2,,l]
  zmin = Prange[1,,l]
  polygon(c(xind, rev(xind)), c(zmax, rev(zmin)),
          col = rgb(0.85,0.85,0.85,0.5), border = NA)
  
  lines(xind,yhat,col='red', ylim=c(lb, ub), lty = 6, type='l',lwd=3)
  lines(xind, y, col='black', ylim=c(lb, ub), lty = 1, type='l',lwd= 3)
  
  line_color = c("#2171B5","#084594")
  #line_style = c("dashed","dotted")
  
  for (l in (1:Alpha_len)){
    
    dim (S_list[l])
    S = S_list[[l]]
    S = S[n,,]
    N = dim(S)[1]
    Nu = dim(S)[2]
    
    if (!is.null(Mean)){
      S = S * matrix(rep(Std,N),nrow = N, byrow=TRUE) + matrix(rep(Mean,N),nrow = N, byrow=TRUE)
    }
    
    PD = matrix(nrow= N, ncol= Nu)
    for (i in 1:N){
      for (j in 1:Nu){
        PD[i,j] = 1 - abs(sum(S[,j] < S[i,j]) -  sum(S[,j] > S[i,j]))/N
      }
    }
    
    rrange = sort(unique(as.vector(PD)))
    rrangeLen = length(rrange)
    dCDF = matrix(nrow= N, ncol= rrangeLen)
    for (i in 1:N){
      for (j in 1:rrangeLen){
        dCDF[i,j] = sum(PD[i,] <= rrange[j])/Nu
      }
    }
    dCDF_original = dCDF
    
    rank_value = vector(length = N)
    for (i in 1:rrangeLen){
      ranked_num =  sum(rank_value!=0)
      if (ranked_num >= 2) {
        cat(i-1,rrangeLen,"\n")
        break
      }
      rank_value[(dCDF[,i] > 0)] = rank(dCDF[, i],ties.method='average')[dCDF[,i] > 0]- ranked_num
      rank_value[round(rank_value) != rank_value] = 0
      ranked_index = (rank_value!=0) 
      dCDF[ranked_index,] = 0
    } 
    
    ED = (N + 1 -  rank_value)/N
    
    # plot median and most extrem 2 curves
    #lines(xind,S[order(ED)[length(ED)],],col='red', ylim=c(lb, ub), type='l',lwd=2)
    lines(xind,S[order(ED)[1],],col=line_color[l], ylim=c(lb, ub), type='l',lwd=4,lty=l+1)
    lines(xind,S[order(ED)[2],],col=line_color[l], ylim=c(lb, ub), type='l',lwd=4,lty=l+1)
  }
  
  
  
  color = vector()
  color[1] = rgb(0.85,0.85,0.85,0.5)
  
  axis(1,at = xind[c(1,13,25,37,47)],labels = xind2[c(1,13,25,37,47)])
  
  legend("top",legend = c("90% predictive region",
                          "true curve","predictive mean curve",
                          "40% quantile curve","99% quantile curve"),
         col= c(NA,"black","red","#2171B5","#084594"),
         lty = c(0,1,6,2,3),
         lwd = c(0,1,1,1,1),
         pch = c(22,NA,NA,NA,NA),
         pt.bg = c(color,NA,NA,NA,NA),
         pt.cex = c(2,NA,NA,NA,NA),
         cex = 0.6)
  title(title_name)
  #mtext(paste(n,Status))
  mtext(Status)
}

ConfidenceIntervalPlot4 =  function(TruePrange, S_list, Prange_all, y_all, n, Nu, Alpha_list, Mean= NULL, Std = NULL,
                                    title_name = "PredictiveRegion", yhat = NULL){
  
  #N = dim(Prange_all)[2]
  
  Alpha_len = length(Alpha_list)
  
  Prange = array(dim=c(dim(Prange_all)[1],dim(Prange_all)[3],dim(Prange_all)[4]))
  
  Prange[,,] = Prange_all[,n,,]
  TruePrange = TruePrange[,n,]
  y = y_all[n,]
  yhat = yhat[n,]
  
  if (!is.null(Mean)){
    Mean2 = matrix(rep(Mean,2),nrow = 2, byrow = TRUE)
    Std2 = matrix(rep(Std,2),nrow = 2, byrow = TRUE)
    for (l in 1:Alpha_len){
      Prange[,,l] = Prange[,,l] * Std2 + Mean2
    }
    y =  y*Std + Mean
    yhat = yhat * Std + Mean
  }
  
  
  
  lb = -2 #min(min(Prange),min(TruePrange))
  ub = 2.2 #max(max(Prange),max(TruePrange))
  
  xind = -1 + (0:(Nu-1))*0.04
  #xind2 = xind
  
  plot(xind, y, xaxt="n", col='orange', ylim=c(lb, ub), type='l',lwd=1,xlab='',ylab='')
  # for (l in rev(1:Alpha_len)){
  #   Alpha = Alpha_list[l]
  #   zmax = Prange[2,,l]
  #   zmin = Prange[1,,l]
  #   polygon(c(xind, rev(xind)), c(zmax, rev(zmin)),
  #           col = rgb(1-Alpha,1-Alpha,1-Alpha,0.7), border = NA)
  # }
  
  l = Alpha_len
  Alpha = Alpha_list[l]
  zmax = Prange[2,,l]
  zmin = Prange[1,,l]
  polygon(c(xind, rev(xind)), c(zmax, rev(zmin)),
          col = rgb(0.80,0.80,0.80,0.7), border = NA)
  
  if((max(y-Prange[2,,l])>0)|(max(Prange[1,,l]-y)>0)){
    Status = "Out of 90% Predictive Region"
  }else{
    Status = "Within 90% Predictive Region"
  }
  
  zmax = TruePrange[2,]
  zmin = TruePrange[1,]
  polygon(c(xind, rev(xind)), c(zmax, rev(zmin)),
          col = rgb(0.82,0.93,0.93,0.5), border = NA) 
  
  lines(xind,yhat,col='red', ylim=c(lb, ub), lty = 6,type='l',lwd=2)
  lines(xind, y, col= 'black', ylim=c(lb, ub), lty = 1,type='l',lwd= 2)
  
  #"#D6604D" "#F4A582" "#FDDBC7"
  #"#084594", "#2171B5", "#4292C6" 
  line_color = c("#2171B5","#084594")
  line_type = c(3,2)
  
  for (l in (1:Alpha_len)){
    S = S_list[[l]]
    S = S[n,,]
    N = dim(S)[1]
    Nu = dim(S)[2]
    
    if (!is.null(Mean)){
      S = S * matrix(rep(Std,N),nrow = N, byrow=TRUE) + matrix(rep(Mean,N),nrow = N, byrow=TRUE)
    }
    
    PD = matrix(nrow= N, ncol= Nu)
    for (i in 1:N){
      for (j in 1:Nu){
        PD[i,j] = 1 - abs(sum(S[,j] < S[i,j]) -  sum(S[,j] > S[i,j]))/N
      }
    }
    
    rrange = sort(unique(as.vector(PD)))
    rrangeLen = length(rrange)
    dCDF = matrix(nrow= N, ncol= rrangeLen)
    for (i in 1:N){
      for (j in 1:rrangeLen){
        dCDF[i,j] = sum(PD[i,] <= rrange[j])/Nu
      }
    }
    dCDF_original = dCDF
    
    rank_value = vector(length = N)
    
    
    for (i in 1:rrangeLen){
      ranked_num =  sum(rank_value!=0)
      if (ranked_num >= 10) {
        cat(i-1,rrangeLen,"\n")
        break
      }
      rank_value[(dCDF[,i] > 0)] = rank(dCDF[, i],ties.method='average')[dCDF[,i] > 0]- ranked_num
      rank_value[round(rank_value) != rank_value] = 0
      ranked_index = (rank_value!=0)
      dCDF[ranked_index,] = 0
    }
    
    ED = (N + 1 -  rank_value)/N
    
    # plot median and most extrem 2 curves
    #lines(xind,S[order(ED)[length(ED)],],col='red', ylim=c(lb, ub), type='l',lwd=2)
    lines(xind,S[order(ED)[1],],col=line_color[l], ylim=c(lb, ub), type='l',lwd=2,lty=l+1)
    lines(xind,S[order(ED)[2],],col=line_color[l], ylim=c(lb, ub), type='l',lwd=2,lty=l+1)
  }
  
  
  transp  = 0.7
  color = vector()

  color[1] = rgb(0.80,0.80,0.80,0.7)
  color[2] = rgb(0.82,0.93,0.93,0.5)
  
  axis(1,at = xind[c(1,26,51)],labels = xind[c(1,26,51)])
  
  legend("top",legend = c("90% predictive region","90% true region",
                          "true curve", "preditive mean curve",
                          "40% quantile curve","90% quantile curve"),
         col= c(rep(NA,1+1),"black","red",line_color),
         lty = c(rep(0,1+1),1,6,2,3),
         lwd = c(rep(0,1+1),1,1,1,1),
         pch = c(rep(22,1+1),NA,NA,NA,NA),
         pt.bg = c(color,NA,NA,NA,NA),
         pt.cex = c(rep(2,1+1),NA,NA,NA,NA),
         cex = 0.6) 
  #title(title_name,font=3)
  #mtext(paste(n,Status)) 
  mtext(Status,font = 0.5)
  
}
 