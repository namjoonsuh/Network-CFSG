####################  preparation  ####################
rm(list=ls())
setwd('C:/Users/Namjoon Suh/Desktop/Stat+PDE/Citation-Network/Codes & Data/Codes')
library('Matrix')      ## For matrix computation on BIC
library('igraph')      ## For drawing graph
source('SynData.R')    ## Function for creating binary network
source('ADMM_Optim.R') ## Function of ADMM algorithm for estimation
source('functions.R')
source('Synthetic Networks.R')
source('GD2.R')
#################### data generation N = 30, K = 3, NNZ = 10 ####################
set.seed(123)
par(mfrow=c(2,3),mar = c(2, 2, 2, 2)) # ??? ??? ??? ??? # 
######## First Scenario ########
set.seed(123); Net1 <- Network1(30,3,9,1)
set.seed(123); Net2 <- Network1(80,4,18,2)
set.seed(123); Net3 <- Network1(120,5,30,3)
######## Second Scenario ########
set.seed(123); Net4 <- Network2(120,10,0,3,20,4)
set.seed(123); Net5 <- Network2(210,0,50,3,18,5)
set.seed(123); Net6 <- Network2(210,10,10,3,18,6)
par(mfrow=c(1,1))

######## Plot Scree ########
source('functions.R')
par(mfrow=c(2,3),mar = c(2, 2, 2, 1))
plotScree(Net1[[1]],1); plotScree(Net2[[1]],2); plotScree(Net3[[1]],3);
plotScree(Net4[[1]],4); plotScree(Net5[[1]],5); plotScree(Net6[[1]],6);
par(mfrow=c(1,1))

################ Scenario 1 ################
gamma = seq(from=0.01,to=0.02,by=0.001);
delta = seq(from=0.02,to=0.03,by=0.001);
Res1 <- Model_Sel(Net1[[1]],gamma,delta) ## gamma : 0.011, delta : 0.023 ##

gamma = seq(from=0.004,to=0.005,by=0.0001);
delta = seq(from=0.009,to=0.01,by=0.0001);
Res2 <- Model_Sel(Net2[[1]],gamma,delta) ## gamma : 0.004, delta : 0.009 ## 

gamma = seq(from=0.002,to=0.003,by=0.0001);
delta = seq(from=0.0046,to=0.0048,by=0.00002);
Res3 <- Model_Sel(Net3[[1]],gamma,delta) ## gamma : 0.003, delta : 0.00475 ##

################ Scenario 2 ################
gamma = seq(from=0.002,to=0.004,by=0.0002);
delta = seq(from=0.015,to=0.02,by=0.0005);
Res4 <- Model_Sel(Net4[[1]],gamma,delta) ## gamma : 0.0032, delta : 0.019 

gamma = seq(from=0.00204,to=0.00244,by=0.00004); 
delta = seq(from=0.045,to=0.05,by=0.0005);
Res5 <- Model_Sel(Net5[[1]],gamma,delta) ## gamma : 0.00212, delta : 0.048

gamma = seq(from=0.002,to=0.004,by=0.0002);
delta = seq(from=0.015,to=0.02,by=0.0005);
Res6 <- Model_Sel(Net6[[1]],gamma,delta) ## gamma : 0.002, delta : 0.016

############################################
source('Synthetic Networks.R')
par(mfrow=c(2,3),mar = c(2, 2, 2, 2))
set.seed(1234); K_means(Net1[[1]],0.011,0.024,3,1); 
set.seed(1234); K_means(Net2[[1]],0.004,0.009,4,2); 
set.seed(1234); K_means(Net3[[1]],0.003,0.00476,5,3);

set.seed(1234); K_means(Net4[[1]],0.0032,0.019,4,4); 
set.seed(1234); K_means(Net5[[1]],0.00212,0.048,5,5); 
set.seed(1234); K_means(Net6[[1]],0.002,0.016,5,6);
par(mfrow=c(1,1))

object_list <- list(Res1,Res2,Res3,Res4,Res5,Res6,Res7,Res8,Res9)
save(object_list,file="myfile.RData")

set.seed(1234); K_means(Net1[[1]],gamma[11],delta[11],2,1); 
set.seed(1234); K_means(Net2[[1]],gamma[1],delta[11],3,2); 
set.seed(1234); K_means(Net3[[1]],gamma[3],delta[5],4,3);
set.seed(1234); K_means(Net4[[1]],0.0032,0.019,4,4); 
set.seed(1234); K_means(Net5[[1]],0.00212,0.048,5,5); 
set.seed(1234); K_means(Net6[[1]],0.002,0.016,5,6);
