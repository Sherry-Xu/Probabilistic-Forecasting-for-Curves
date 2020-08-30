corenum = 24
N_rep = 400
DataType = "Sim"
runtime = "20200830"

threshold = 0.999; Alpha_list = c(0.60,0.10); #0.1
Ncv = 200; Jcv = 6
Ns = 200; Nu=51; 
Y = NULL;

p = 1; N = 100; D = 4; sigmaE = 0.25; dist = 'norm'  #"t","exp","norm"
#Kchi2 = 2000 

Kchi2 = 1500

HO = TRUE;
sparsity = FALSE #"diag","log", FALSE
  
ifplot = FALSE
  