rm(list=ls())
setwd('C:/Users/namjo/Desktop/Research Files/Network Paper/Codes & Data/Codes')
library('Matrix')      ## For matrix computation on BIC
library('igraph')      ## For drawing graph
library('extrafont')

source('Synthetic Networks.R')
source('SynData.R')    ## Function for creating binary network
source('ADMM_Optim.R') ## Function of ADMM algorithm for estimation
source('ADMM2_Optim.R') ## Function of ADMM algorithm for estimation
source('ADMM3_Optim.R')
source('functions.R')
source('Synthetic Networks.R')

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

################### ADMM 1 
################### For fixed sparsity ###################
######################### K = 4 #########################
gamma = seq(from=0.00455,to=0.00460,by=0.00001)
delta = seq(from=0.000,to=0.015,by=0.0005)

Case1 <- Model_Sel(Net1[[1]],gamma,delta);
Res1 <- ADMM(Net1[[1]],0.00455,0.011)
qr(Res1[[3]])$rank

### Frobenius Norm difference : 0.2156143
sum((Res1[[3]]-Normalize(Net1[[4]]))^2)

######################### K = 5 #########################
gamma = seq(from=0.00455,to=0.00460,by=0.00001)
delta = seq(from=0.008,to=0.015,by=0.001)

Case2 <- Model_Sel(Net2[[1]],gamma,delta);
Res2 <- ADMM(Net2[[1]],0.00458,0.009)
qr(Res2[[3]])$rank;

### Frobenius Norm difference : 0.3254004
sum((Res2[[3]]-Normalize(Net2[[4]]))^2)


######################### K = 6 #########################
gamma = seq(from=0.00465,to=0.00470,by=0.00001)
delta = seq(from=0.008,to=0.014,by=0.001)

Case3 <- Model_Sel(Net3[[1]],gamma,delta);
Res3 <- ADMM(Net3[[1]],0.00468,0.009)
qr(Res3[[3]])$rank

### Frobenius Norm difference : 0.4484387
sum((Res3[[3]]-Normalize(Net3[[4]]))^2)


######################### K = 7 #########################
gamma = seq(from=0.0047,to=0.00475,by=0.00001)
delta = seq(from=0.007,to=0.008,by=0.0002)

Case4 <- Model_Sel(Net4[[1]],gamma,delta);
Res4 <- ADMM(Net4[[1]],0.00475,0.007)
qr(Res4[[3]])$rank

### Frobenius Norm difference : 0.600115
sum((Res4[[3]]-Normalize(Net4[[4]]))^2)

Diff <- c(0.2156143,0.3254004,0.4484387,0.600115)
plot(1:4,c(93.46956,93.95154,94.47971,94.80496))

################### ADMM 2 ###################
################### For fixed sparsity ###################
######################### K = 4 #########################
gamma = seq(from=0.00455,to=0.00460,by=0.00001)
delta = seq(from=0.011,to=0.0115,by=0.0001)

Case1 <- Model_Sel2(Net1[[1]],gamma,delta);
Res1 <- ADMM2(Net1[[1]],0.00455,0.011)
qr(Res1[[3]])$rank

### Frobenius Norm difference : 113.7751
sum((Normalize(Res1[[4]])-Normalize(Net1[[5]]))^2)

######################### K = 5 #########################
gamma = seq(from=0.00455,to=0.00460,by=0.00001)
delta = seq(from=0.008,to=0.015,by=0.001)

Case2 <- Model_Sel2(Net2[[1]],gamma,delta);
Res2 <- ADMM2(Net2[[1]],0.00458,0.009)
qr(Res2[[3]])$rank;

### Frobenius Norm difference : 0.3245855
sum((Normalize(Res2[[3]])-Normalize(Net2[[5]]))^2)


######################### K = 6 #########################
gamma = seq(from=0.00465,to=0.00470,by=0.00001)
delta = seq(from=0.008,to=0.014,by=0.001)

Case3 <- Model_Sel2(Net3[[1]],gamma,delta);
Res3 <- ADMM2(Net3[[1]],0.00468,0.009)
qr(Res3[[3]])$rank

### Frobenius Norm difference : 0.4479127
sum((Normalize(Res3[[3]])-Normalize(Net3[[4]]))^2)


######################### K = 7 #########################
gamma = seq(from=0.0047,to=0.00475,by=0.00001)
delta = seq(from=0.007,to=0.008,by=0.0002)

Case4 <- Model_Sel(Net4[[1]],gamma,delta);
Res4 <- ADMM2(Net4[[1]],0.00475,0.007)
qr(Res4[[3]])$rank

### Frobenius Norm difference : 0.600115
sum((Normalize(Res4[[3]])-Normalize(Net4[[4]]))^2)



################### ADMM 1 

################### For fixed rank of L ###################
######################### |S| = 18 #########################
gamma = seq(from=0.00455,to=0.00460,by=0.00001)
delta = seq(from=0.010,to=0.013,by=0.0005)

Case11 <- Model_Sel(Net11[[1]],gamma,delta);
Res11 <- ADMM(Net11[[1]],0.00457,0.0105)
non_zero_S <- 0
for(i in 1:100){
  for(j in 1:100){
    if(i<j){
      if(Res11[[4]][i,j] > 0){
        non_zero_S = non_zero_S + 1
      }
    }
  }
}
non_zero_S

### Frobenius Norm difference : 14110.3
sum((Normalize(Res11[[4]])-Normalize(Net11[[5]]))^2)
sum((Res11[[4]]-Net11[[5]])^2)


######################### |S| = 30 #########################
gamma = seq(from=0.00455,to=0.00460,by=0.00001)
delta = seq(from=0.010,to=0.013,by=0.0005)

Case22 <- Model_Sel(Net22[[1]],gamma,delta);
Res22 <- ADMM(Net22[[1]],0.00460,0.0105)
non_zero_S <- 0
for(i in 1:100){
  for(j in 1:100){
    if(i<j){
      if(Res22[[4]][i,j] > 0){
        non_zero_S = non_zero_S + 1
      }
    }
  }
}
non_zero_S

### Frobenius Norm difference : 22058.6
sum((Normalize(Res22[[4]])-Normalize(Net22[[5]]))^2)
sum((Res22[[4]]-Net22[[5]])^2)

######################### |S| = 42 #########################
gamma = seq(from=0.00455,to=0.00460,by=0.00001)
delta = seq(from=0.010,to=0.013,by=0.0005)

Case33 <- Model_Sel(Net33[[1]],gamma,delta);
Res33 <- ADMM(Net33[[1]],0.00459,0.0115)
non_zero_S <- 0
for(i in 1:100){
  for(j in 1:100){
    if(i<j){
      if(Res33[[4]][i,j] > 0){
        non_zero_S = non_zero_S + 1
      }
    }
  }
}
non_zero_S

### Frobenius Norm difference : 32506.73
sum((Normalize(Res33[[4]])-Normalize(Net33[[5]]))^2)
sum((Res33[[4]]-Net33[[5]])^2)

######################### |S| = 54 #########################
gamma = seq(from=0.00455,to=0.00460,by=0.00001)
delta = seq(from=0.010,to=0.013,by=0.0005)

Case44 <- Model_Sel(Net44[[1]],gamma,delta);
Res44 <- ADMM(Net44[[1]],0.0046,0.0105)
non_zero_S <- 0
for(i in 1:100){
  for(j in 1:100){
    if(i<j){
      if(Res44[[4]][i,j] > 0){
        non_zero_S = non_zero_S + 1
      }
    }
  }
}
non_zero_S

### Frobenius Norm difference : 42821.23
sum((Res44[[4]]-Net44[[5]])^2)
plot(1:4,c(14110.3,22058.6,32506.73,42821.23))
 
################### ADMM 2 
################### For fixed rank of L ###################
######################### |S| = 18 #########################
gamma = seq(from=0.00455,to=0.00460,by=0.00001)
delta = seq(from=0.010,to=0.013,by=0.0005)

Case11 <- Model_Sel2(Net11[[1]],gamma,delta);
Res11 <- ADMM2(Net11[[1]],0.00457,0.0105)
qr(Res11[[3]])$rank
non_zero_S <- 0
for(i in 1:100){
  for(j in 1:100){
    if(i<j){
      if(Res11[[4]][i,j] > 0){
        non_zero_S = non_zero_S + 1
      }
    }
  }
}
non_zero_S

### Frobenius Norm difference : 14110.3
sum((Res11[[4]]-Net11[[5]])^2)

######################### |S| = 30 #########################
gamma = seq(from=0.00455,to=0.00460,by=0.00001)
delta = seq(from=0.010,to=0.013,by=0.0005)

Case22 <- Model_Sel2(Net22[[1]],gamma,delta);
Res22 <- ADMM2(Net22[[1]],0.00460,0.0105)
qr(Res22[[3]])$rank
non_zero_S <- 0
for(i in 1:100){
  for(j in 1:100){
    if(i<j){
      if(Res22[[4]][i,j] > 0){
        non_zero_S = non_zero_S + 1
      }
    }
  }
}
non_zero_S

### Frobenius Norm difference : 22058.59
sum((Res22[[4]]-Net22[[5]])^2)

######################### |S| = 42 #########################
gamma = seq(from=0.00455,to=0.00460,by=0.00001)
delta = seq(from=0.010,to=0.013,by=0.0005)

Case33 <- Model_Sel2(Net33[[1]],gamma,delta);
Res33 <- ADMM2(Net33[[1]],0.00459,0.0115)
qr(Res33[[3]])$rank
non_zero_S <- 0
for(i in 1:100){
  for(j in 1:100){
    if(i<j){
      if(Res33[[4]][i,j] > 0){
        non_zero_S = non_zero_S + 1
      }
    }
  }
}
non_zero_S

### Frobenius Norm difference : 32506.67
sum((Res33[[4]]-Net33[[5]])^2)

######################### |S| = 54 #########################
gamma = seq(from=0.00455,to=0.00460,by=0.00001)
delta = seq(from=0.010,to=0.013,by=0.0005)

Case44 <- Model_Sel(Net44[[1]],gamma,delta);
Res44 <- ADMM2(Net44[[1]],0.0046,0.0105)
non_zero_S <- 0
for(i in 1:100){
  for(j in 1:100){
    if(i<j){
      if(Res44[[4]][i,j] > 0){
        non_zero_S = non_zero_S + 1
      }
    }
  }
}
non_zero_S

### Frobenius Norm difference : 42821.26
sum((Res44[[4]]-Net44[[5]])^2)


################### Draw plot ################### 
par(mfrow=c(1,2))
Diff1 <- t(t(c(0.2156143,0.3254004,0.4484387,0.600115)))
Diff2 <- t(t(c(0.2154249,0.3245855,0.4479127,0.600115)))
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
Diff1 <- t(t(c(14110.3,22058.6,32506.73,42821.23)))
Diff2 <- t(t(c(14110.3,22058.59,32506.67,42821.26)))
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

