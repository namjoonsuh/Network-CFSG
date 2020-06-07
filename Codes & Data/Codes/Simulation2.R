rm(list=ls())
setwd('C:/Users/namjo/Desktop/Research Files/Network Paper/Codes & Data/Codes')
library('Matrix')      ## For matrix computation on BIC
library('igraph')      ## For drawing graph
library('extrafont')

source('Synthetic Networks.R')
source('SynData.R')     ## Function for creating binary network
source('ADMM_Optim.R')  ## Function of ADMM algorithm for estimation  # delta = 10^-4, r = 15
source('ADMM2_Optim.R') ## Function of ADMM algorithm for estimation

Normalize <- function(X){
  N = nrow(X); Norm = matrix(0,nrow=N,ncol=N);
  for(i in 1:N){
    if(sum(X[,i])==0){
      next;
    }
    Norm[,i] = X[,i]/sqrt(sum(X[,i]^2))
  }
  return(Norm)
}
Non_Zero_S <- function(X){
  N = nrow(X); count <- 0;
  for(i in 1:N){
    for(j in 1:N){
      if(i<j && X[i,j]!=0){
        count <- count + 1;
      }
    }
  }
  print(count)
}

############### Create Network for simulation ###############
font_import()
fonts()
loadfonts(device="win")

par(mfrow=c(2,4))
set.seed(123); Net1 <- Network1(100,4,50,1);
set.seed(123); Net2 <- Network1(100,5,50,1);
set.seed(123); Net3 <- Network1(100,6,50,1);
set.seed(123); Net4 <- Network1(100,7,50,1);

set.seed(123); Net11 <- Network1(100,4,18,1);
set.seed(123); Net22 <- Network1(100,4,30,1);
set.seed(123); Net33 <- Network1(100,4,42,1);
set.seed(123); Net44 <- Network1(100,4,54,1);

##########################################################
###################### Experiments #######################

Res1 <- ADMM(Net1[[1]],0.004,0.0085)
qr(Res1[[3]])$rank; Non_Zero_S(Res1[[4]]);
Res2 <- ADMM(Net2[[1]],0.004,0.008)
qr(Res2[[3]])$rank; Non_Zero_S(Res2[[4]]);
Res3 <- ADMM(Net3[[1]],0.004,0.007)
qr(Res3[[3]])$rank; Non_Zero_S(Res3[[4]]);
Res4 <- ADMM(Net4[[1]],0.0045,0.007)
qr(Res4[[3]])$rank; Non_Zero_S(Res4[[4]]);


Res11 <- ADMM(Net11[[1]],0.0038,0.008)
qr(Res11[[3]])$rank; Non_Zero_S(Res11[[4]]);
Res22 <- ADMM(Net22[[1]],0.0037,0.0084)
qr(Res22[[3]])$rank; Non_Zero_S(Res22[[4]]);
Res33 <- ADMM(Net33[[1]],0.0037,0.0091)
qr(Res33[[3]])$rank; Non_Zero_S(Res33[[4]]);
Res44 <- ADMM(Net44[[1]],0.0035,0.0085);
qr(Res44[[3]])$rank; Non_Zero_S(Res44[[4]]);

################### ADMM 1 ###############
################### For fixed sparsity ###################
######################### K = 4 #########################
gamma = seq(from=0.002,to=0.004,by=0.001)
delta = seq(from=0.007,to=0.010,by=0.0005)

Case1 <- Model_Sel(Net1[[1]],gamma,delta);

Res1 <- ADMM(Net1[[1]],0.004,0.007)
qr(Res1[[3]])$rank; Non_Zero_S(Res1[[4]]);

### Frobenius Norm difference : 0.2451939
sum((Normalize(Res1[[3]])-Normalize(Net1[[4]]))^2)

######################### K = 5 #########################
gamma = seq(from=0.004,to=0.0045,by=0.0002)
delta = seq(from=0.003,to=0.008,by=0.0005)

Case2 <- Model_Sel(Net2[[1]],gamma,delta);
Res2 <- ADMM(Net2[[1]],0.0042,0.0055)
qr(Res2[[3]])$rank; Non_Zero_S(Res2[[4]]);

### Frobenius Norm difference : 0.3718837
sum((Normalize(Res2[[3]])-Normalize(Net2[[4]]))^2)

######################### K = 6 #########################
gamma = seq(from=0.004,to=0.0045,by=0.0005)
delta = seq(from=0.007,to=0.008,by=0.0005)

Case3 <- Model_Sel(Net3[[1]],gamma,delta);
Res3 <- ADMM(Net3[[1]],0.004,0.007)
qr(Res3[[3]])$rank; Non_Zero_S(Res3[[4]]);

### Frobenius Norm difference : 0.5074376
sum((Normalize(Res3[[3]])-Normalize(Net3[[4]]))^2)

######################### K = 7 #########################
gamma = seq(from=0.0045,to=0.00455,by=0.00005)
delta = seq(from=0.007,to=0.008,by=0.0003)

Case4 <- Model_Sel(Net4[[1]],gamma,delta);
Res4 <- ADMM(Net4[[1]],0.0045,0.007)
qr(Res4[[3]])$rank; Non_Zero_S(Res4[[4]]);

### Frobenius Norm difference : 0.7385205
sum((Normalize(Res4[[3]])-Normalize(Net4[[4]]))^2)

Diff <- c(0.2451939,0.3728151,0.5621203,0.7339069)

################### ADMM 2 ###################
################### For fixed sparsity ###################
######################### K = 4 #########################
gamma = seq(from=0.002,to=0.004,by=0.001)
delta = seq(from=0.007,to=0.010,by=0.0005)

Case1 <- Model_Sel2(Net1[[1]],gamma,delta);

Res1 <- ADMM2(Net1[[1]],gamma[3],delta[4])
qr(Res1[[3]])$rank; Non_Zero_S(Res1[[4]]);

### Frobenius Norm difference : 0.2446366
sum((Normalize(Res1[[3]])-Normalize(Net1[[4]]))^2)

######################### K = 5 #########################
gamma = seq(from=0.004,to=0.0045,by=0.0005)
delta = seq(from=0.008,to=0.015,by=0.003)

Case2 <- Model_Sel2(Net2[[1]],gamma,delta);
Res2 <- ADMM2(Net2[[1]],gamma[1],delta[1])
qr(Res2[[3]])$rank; Non_Zero_S(Res2[[4]]);

### Frobenius Norm difference : 0.3719354
sum((Normalize(Res2[[3]])-Normalize(Net2[[4]]))^2)

######################### K = 6 #########################
gamma = seq(from=0.004,to=0.0045,by=0.0005)
delta = seq(from=0.0055,to=0.0065,by=0.0005)

Case3 <- Model_Sel2(Net3[[1]],gamma,delta);
Res3 <- ADMM2(Net3[[1]],gamma[1],delta[1])
qr(Res3[[3]])$rank; Non_Zero_S(Res3[[4]]);

### Frobenius Norm difference : 0.515301
sum((Normalize(Res3[[3]])-Normalize(Net3[[4]]))^2)

######################### K = 7 #########################
gamma = seq(from=0.0045,to=0.00455,by=0.00005)
delta = seq(from=0.0055,to=0.006,by=0.0003)

Case4 <- Model_Sel2(Net4[[1]],gamma,delta);
Res4 <- ADMM2(Net4[[1]],gamma[1],delta[2])
qr(Res4[[3]])$rank; Non_Zero_S(Res4[[4]]);

### Frobenius Norm difference : 0.7367166
sum((Normalize(Res4[[3]])-Normalize(Net4[[4]]))^2)

Diff <- c(0.2086144,0.3728151,0.5621077,0.7387447)

################### ADMM 1 






################### ADMM 1 ###############
################### For fixed rank of L ###################
######################### |S| = 18 #########################
gamma = seq(from=0.003,to=0.004,by=0.0002)
delta = seq(from=0.007,to=0.01,by=0.0005)

Case11 <- Model_Sel(Net11[[1]],gamma,delta);
Res11 <- ADMM(Net11[[1]],0.0038,0.008)
qr(Res11[[3]])$rank; Non_Zero_S(Res11[[4]]);

### Frobenius Norm difference : 14105.93
sum((Res11[[4]]-Net11[[5]])^2)

######################### |S| = 30 #########################
gamma = seq(from=0.0035,to=0.004,by=0.0002)
delta = seq(from=0.007,to=0.01,by=0.0007)

Case22 <- Model_Sel(Net22[[1]],gamma,delta);
Res22 <- ADMM(Net22[[1]],0.0037,0.0084)
qr(Res22[[3]])$rank; Non_Zero_S(Res22[[4]]);

### Frobenius Norm difference : 22051.4
sum((Normalize(Res22[[4]])-Normalize(Net22[[5]]))^2)
sum((Res22[[4]]-Net22[[5]])^2)

######################### |S| = 42 #########################
gamma = seq(from=0.0035,to=0.004,by=0.0002)
delta = seq(from=0.007,to=0.01,by=0.0007)

Case33 <- Model_Sel(Net33[[1]],gamma,delta);
Res33 <- ADMM(Net33[[1]],0.0037,0.0091)
qr(Res33[[3]])$rank; Non_Zero_S(Res33[[4]]);

### Frobenius Norm difference : 32496.14
sum((Normalize(Res33[[4]])-Normalize(Net33[[5]]))^2)
sum((Res33[[4]]-Net33[[5]])^2)

######################### |S| = 54 #########################
gamma = seq(from=0.0035,to=0.004,by=0.0002)
delta = seq(from=0.007,to=0.01,by=0.0005)

Case44 <- Model_Sel(Net44[[1]],gamma,delta);
Res44 <- ADMM(Net44[[1]],0.0035,0.0085);
qr(Res44[[3]])$rank; Non_Zero_S(Res44[[4]]);

### Frobenius Norm difference : 42801.32
sum((Res44[[4]]-Net44[[5]])^2)
plot(1:4,c(14110.3,22058.6,32506.73,42821.23))
 
################### ADMM 2 


################### ADMM 2 ###############
################### For fixed rank of L ###################
######################### |S| = 18 ########################
gamma = seq(from=0.003,to=0.004,by=0.0002)
delta = seq(from=0.006,to=0.01,by=0.0005)

Case11 <- Model_Sel2(Net11[[1]],gamma,delta);
Res11 <- ADMM2(Net11[[1]],gamma[5],delta[5])
qr(Res11[[3]])$rank; Non_Zero_S(Res11[[4]]);

### Frobenius Norm difference : 14106.53
sum((Normalize(Res11[[4]])-Normalize(Net11[[5]]))^2)
sum((Res11[[4]]-Net11[[5]])^2)

######################### |S| = 30 #########################
gamma = seq(from=0.0035,to=0.004,by=0.0002)
delta = seq(from=0.007,to=0.01,by=0.0007)

Case22 <- Model_Sel2(Net22[[1]],gamma,delta);
Res22 <- ADMM2(Net22[[1]],gamma[2],delta[3])
qr(Res22[[3]])$rank; Non_Zero_S(Res22[[4]]);

### Frobenius Norm difference : 22051.44
sum((Normalize(Res22[[4]])-Normalize(Net22[[5]]))^2)
sum((Res22[[4]]-Net22[[5]])^2)

######################### |S| = 42 #########################
gamma = seq(from=0.0035,to=0.004,by=0.0002)
delta = seq(from=0.007,to=0.01,by=0.0007)

Case33 <- Model_Sel2(Net33[[1]],gamma,delta);
Res33 <- ADMM2(Net33[[1]],gamma[2],delta[4])
qr(Res33[[3]])$rank; Non_Zero_S(Res33[[4]]);

### Frobenius Norm difference : 32496.54
sum((Normalize(Res33[[4]])-Normalize(Net33[[5]]))^2)
sum((Res33[[4]]-Net33[[5]])^2)

######################### |S| = 54 #########################
gamma = seq(from=0.0035,to=0.004,by=0.0002)
delta = seq(from=0.007,to=0.01,by=0.0005)

Case44 <- Model_Sel2(Net44[[1]],gamma,delta);
Res44 <- ADMM2(Net44[[1]],gamma[2],delta[4]);
qr(Res44[[3]])$rank; Non_Zero_S(Res44[[4]]);

### Frobenius Norm difference : 42807.57
sum((Res44[[4]]-Net44[[5]])^2)
plot(1:4,c(14110.3,22058.6,32506.73,42821.23))






######################     END     #######################
##########################################################

################### Draw plot ################### 
par(mfrow=c(1,2))
Diff1 <- t(t(c(0.2451939,0.3718837,0.5074376,0.7385205)))
Diff2 <- t(t(c(0.2446366,0.3719354,0.515301,0.7367166)))
value <- rbind(Diff1,Diff2)
Algorithm <- t(t(c("ADMM 1","ADMM 1","ADMM 1","ADMM 1","ADMM 2","ADMM 2","ADMM 2","ADMM 2")))

library(dplyr)
library(ggplot2)
rank <- t(t(seq(4, 7, by=1)))
rank <- rbind(rank,rank)

df <- tbl_df(data.frame(value,Algorithm))

library(ggplot2)
p1<-ggplot(df, aes(x = rank, y = value, shape = Algorithm)) + 
  geom_point(size=3) +
  geom_smooth(method="lm",se = FALSE) +
  ylab(TeX('Squared Frobenius norm of $\\hat{\\Delta}_{N}^{L}$')) 


#### Plot
Diff1 <- t(t(c(14105.93,22051.4,32496.14,42801.32)))
Diff2 <- t(t(c(14106.53,22051.44,32496.54,42807.57)))
value <- rbind(Diff1,Diff2)
Algorithm <- t(t(c("ADMM 1","ADMM 1","ADMM 1","ADMM 1","ADMM 2","ADMM 2","ADMM 2","ADMM 2")))

library(dplyr)
library(ggplot2)
sparsity <- t(t(c(18,30,42,54)))
sparsity <- rbind(sparsity,sparsity)

df <- tbl_df(data.frame(value,Algorithm))
loadfonts(device="win")

library(latex2exp)
library(ggplot2)
p2 <- ggplot(df, aes(x = sparsity, y = value, shape = Algorithm)) + 
  geom_point(size=3) +
  geom_smooth(method="lm",se = FALSE) +
  ylab(TeX('Squared Frobenius norm of $\\hat{\\Delta}^{S}$')) 


library(cowplot)
plot_grid(p1,p2,labels="AUTO")
#### ratio -> 4.5:11 

