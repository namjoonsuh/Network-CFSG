####################  preparation  ####################
rm(list=ls())
setwd('C:/Users/namjo/Documents/GitHub/Citation-Network/Codes & Data/Codes')
library('Matrix')      ## For matrix computation on BIC
library('igraph')      ## For drawing graph
source('SynData.R')    ## Function for creating binary network
source('ADMM_Optim.R') ## Function of ADMM algorithm for estimation
source('functions.R')
source('Synthetic Networks.R')
#################### data generation N = 30, K = 3, NNZ = 10 ####################
set.seed(123)
par(mfrow=c(2,3),mar = c(1, 1.5, 1.5, 1.5)) # ??? ??? ??? ??? # 
######## First Scenario ########
set.seed(123); Net1 <- Network1(30,3,9,1)
set.seed(123); Net2 <- Network1(80,4,18,2)
set.seed(123); Net3 <- Network1(120,5,30,3)
######## Second Scenario ########
set.seed(123); Net4 <- Network2(120,10,0,3,20,4)
set.seed(123); Net5 <- Network2(210,0,50,3,18,5)
set.seed(123); Net6 <- Network2(210,10,10,3,18,6)
par(mfrow=c(1,1))

######## Third Scenario ########
set.seed(123); Net7 <- Network3(120,10,0,3,20,7)
set.seed(123); Net8 <- Network3(210,0,50,3,18,8)
set.seed(123); Net9 <- Network3(210,10,10,3,18,9);
par(mfrow=c(1,1))

######## Plot Scree ########
par(mfrow=c(3,3),mar = c(2, 2, 2, 1))
plotScree(Net1,1); plotScree(Net2,2); plotScree(Net3,3);
plotScree(Net4,4); plotScree(Net5,5); plotScree(Net6,6);
plotScree(Net7,7); plotScree(Net8,8); plotScree(Net9,9);
par(mfrow=c(1,1))

################ Scenario 1 ################
gamma = seq(from=0.01,to=0.02,by=0.001);
delta = seq(from=0.02,to=0.03,by=0.001);
Res1 <- Model_Sel(Net1,gamma,delta) ## gamma : 0.011, delta : 0.023 

gamma = seq(from=0.004,to=0.005,by=0.0001);
delta = seq(from=0.009,to=0.01,by=0.0001);
Res2 <- Model_Sel(Net2,gamma,delta) ## gamma : 0.004, delta : 0.009 ## 

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

################ Scenario 3 ################
gamma = seq(from=0.002,to=0.004,by=0.0002);
delta = seq(from=0.015,to=0.02,by=0.0005);
Res7 <- Model_Sel(Net7[[1]],gamma,delta) ## gamma : 0.0032, delta : 0.019 
Eval7 <- Eval(Net7,Res7,gamma,delta)

gamma = seq(from=0.0019,to=0.002,by=0.00001);
delta = seq(from=0.03,to=0.035,by=0.0005);
Res8 <- Model_Sel(Net8[[1]],gamma,delta) ## gamma : 0.00212, delta : 0.048

gamma = seq(from=0.001,to=0.002,by=0.0001);
delta = seq(from=0.022,to=0.024,by=0.0002);
Res9 <- Model_Sel(Net9[[1]],gamma,delta) ## gamma : 0.002, delta : 0.016
Eval9 <- Eval(Net9,Res9,gamma,delta)
  
############################################
source('Synthetic Networks.R')
par(mfrow=c(3,3),mar = c(2, 2, 2, 1))
set.seed(1234); K_means(Net1[[1]],0.011,0.024,3,1); 
set.seed(1234); K_means(Net2[[1]],0.004,0.009,4,2); 
set.seed(1234); K_means(Net3[[1]],0.003,0.00476,5,3);

set.seed(1234); K_means(Net4[[1]],0.0032,0.019,3,4); 
set.seed(1234); K_means(Net5[[1]],0.00212,0.048,5,5); 
set.seed(1234); K_means(Net6[[1]],0.002,0.016,5,6);

set.seed(1234); K_means(Net7[[1]],0.0032,0.019,4,7); 
set.seed(1234); K_means(Net8[[1]],0.00192,0.0315,6,8); 
set.seed(12345); K_means(Net9[[1]],0.0019,0.024,6,9);
par(mfrow=c(1,1))

object_list <- list(Res1,Res2,Res3,Res4,Res5,Res6,Res7,Res8,Res9)
save(object_list,file="myfile.RData")

(sum(Net1)/2)*0.03
(sum(Net2)/2)*0.01
(sum(Net3)/2)*0.01

(sum(Net4[[1]])/2)*0.001
(sum(Net5[[1]])/2)*0.001
(sum(Net6[[1]])/2)*0.001

(sum(Net7[[1]])/2)*0.001
(sum(Net8[[1]])/2)*0.001
(sum(Net9[[1]])/2)*0.001



library(extrafont)
loadfonts(device = "win")
par(family = "LM Roman 10")
x <- seq(1, 10, 1)
y <- seq(1, 10, 1)
plot(y ~ x, main="This plot uses LaTeX font!", ylab = expression(alpha))


library(cluster)
data(votes.repub)
agn1 <- agnes(votes.repub, metric = "manhattan", stand = TRUE)
plot(agn1, which.plots=2)
