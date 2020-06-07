rm(list=ls())
setwd('C:/Users/namjo/Desktop/Research Files/Network Paper/Codes & Data/Codes')
library(latex2exp)
source('ADMM_Optim.R') # delta = 10^-4, r = 15
source('functions.R')
source('Synthetic Networks.R')

CrtData <- function(alpha,L,S,N){
  # Create empty matrix for storing adjacency matrix X and Probability for 
  X <- matrix(0,N,N); 
  P <- matrix(0,N,N);
  P <- exp(alpha+L+S)/(1+exp(alpha+L+S));
  
  # Create random graph according to (a, L, S)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      if(P[i,j]>runif(1,0,1))
        X[i,j]<-1
    }
  }
  X <- X + t(X)
  result <- list(X,P)
  return(result)
}
CV <- function(X,gamma,delta,K){
  rate <- 0;
  for(j in 1:K){
    X_fit = X; count = 0;
    n = nrow(X); I_1 = sample(1:n, floor(n/2)); I_2 = setdiff(1:n,I_1);
    M = expand.grid(I_2,I_2);
    for(i in 1:length(I_2)^2){ X_fit[M[i,1],M[i,2]]=0 }
    res = ADMM(X_fit,gamma,delta)
    alpha = res[[1]]; L = res[[3]]; S = res[[4]];
    X_new = CrtData(alpha,L,S,n)[[1]];
    
    for(l in 1:length(I_2)^2){ 
      if(X[M[l,1],M[l,2]]!=X_new[M[l,1],M[l,2]])
        count <- count + 1;
    }
    rate <- rate + (count/length(I_2)^2)/K;
  }
  return(rate)
}

######## Draw figure 1 in paper 
######## Generate network and draw figures of networks for 6 cases ########
set.seed(123)
par(mfrow=c(2,3),mar = c(2, 2, 2, 2)) 
######## First Scenario ########
set.seed(123); Net1 <- Network1(30,3,9,1)
set.seed(123); Net2 <- Network1(80,4,18,2)
set.seed(123); Net3 <- Network1(120,5,30,3)

######## Second Scenario ########
set.seed(123); Net4 <- Network2(120,10,0,3,20,4)
set.seed(123); Net5 <- Network2(210,0,50,3,18,5)
set.seed(123); Net6 <- Network2(210,10,10,3,18,6)
par(mfrow=c(1,1))

######## Draw figure 2 in paper (scree Plot of each case) ########
par(mfrow=c(2,3),mar = c(2, 2, 2, 2)) 
plotScree(Net1[[1]],1); plotScree(Net2[[1]],2); plotScree(Net3[[1]],3); 
plotScree(Net4[[1]],4); plotScree(Net5[[1]],5); plotScree(Net6[[1]],6); 
par(mfrow=c(1,1))

######################### Scenario 1 #########################
######################### Case 1 #########################
gamma = seq(from=0.01,to=0.02,by=0.002);
delta = seq(from=0.02,to=0.03,by=0.002);

Case1 <- Model_Sel(Net1[[1]],gamma,delta);

count <- 1;
MisCl_rate1 <- matrix(rep(0,length(gamma)*length(delta)),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate1[i,j] <- CV(Net1[[1]],gamma[i],delta[j],10)    
    print(count)
    count <- count + 1;
  }
}

CV_ind <- which(MisCl_rate1==min(MisCl_rate1[Case1[[4]]==3]),arr.ind=TRUE)

Eval_func(Net1,gamma[4],delta[6],1) # AIC
Eval_func(Net1,gamma[4],delta[6],1) # BIC
Eval_func(Net1,gamma[5],delta[2],1) # NHCV

######################### Case 2 #########################
gamma = seq(from=0.004,to=0.005,by=0.0002)
delta = seq(from=0.0000,to=0.0002,by=0.00005)

Case2 <- Model_Sel(Net2[[1]],gamma,delta);

count <- 1;
MisCl_rate2 <- matrix(rep(0,length(gamma)^2),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate2[i,j] <- CV(Net2[[1]],gamma[i],delta[j],10)    
    print(count)
    count <- count + 1;
  }
}
CV_ind <- which(MisCl_rate2==min(MisCl_rate2[Case2[[4]]==4]),arr.ind=TRUE)

Eval_func(Net2,gamma[4],delta[6],2) # AIC
Eval_func(Net2,gamma[4],delta[6],2) # BIC
Eval_func(Net2,gamma[4],delta[2],2) # NHCV

######################### Case 3 #########################
gamma = seq(from=0.003,to=0.004,by=0.0002)
delta = seq(from=0.004,to=0.0044,by=0.00004)

Case3 <- Model_Sel(Net3[[1]],gamma,delta);

count <- 1;
MisCl_rate3 <- matrix(rep(0,length(gamma)^2),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate3[i,j] <- CV(Net3[[1]],gamma[i],delta[j],10)    
    print(count)
    count <- count + 1;
  }
}

CV_ind <- which(MisCl_rate3==min(MisCl_rate3[Case3[[4]]==5]),arr.ind=TRUE)

Eval_func(Net3,gamma[1],delta[4],3) # AIC
Eval_func(Net3,gamma[1],delta[4],3) # BIC
Eval_func(Net3,gamma[5],delta[7],3) # NHCV

######################### Scenario 2 #########################
######################### Case 4 #########################
gamma = seq(from=0.003,to=0.004,by=0.0002);
delta = seq(from=0.055,to=0.065,by=0.002);
Case4 <- Model_Sel(Net4[[1]],gamma,delta);

count <- 1;
MisCl_rate4 <- matrix(rep(0,length(gamma)^2),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate4[i,j] <- CV(Net4[[1]],gamma[i],delta[j],10)    
    print(count)
    count <- count + 1;
  }
}

CV_ind <- which(MisCl_rate4==min(MisCl_rate4[Case4[[4]]==3]),arr.ind=TRUE)

Eval_func(Net4,gamma[4],delta[5],4) # AIC
Eval_func(Net4,gamma[4],delta[5],4) # BIC
Eval_func(Net4,gamma[3],delta[1],4) # NHCV

K_means(Net4[[1]], gamma[3],delta[1], 4, 4)
sort(Net4[[7]]) # Nodes with 3 overlapped membership

######################### Case 5 #########################
gamma = seq(from=0.002,to=0.0022,by=0.00005);
delta = seq(from=0.045,to=0.055,by=0.0025);
Case5 <- Model_Sel(Net5[[1]],gamma,delta);

## Network CV
count <- 1;
MisCl_rate5 <- matrix(rep(0,length(gamma)^2),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate5[i,j] <- CV(Net5[[1]],gamma[i],delta[j],10)    
    print(count)
    count <- count + 1;
  }
}

CV_ind <- which(MisCl_rate5==min(MisCl_rate5[Case5[[4]]==3]),arr.ind=TRUE)

Eval_func(Net5,gamma[5],delta[1],5) # AIC
Eval_func(Net5,gamma[5],delta[2],5) # BIC
Eval_func(Net5,gamma[4],delta[2],5) # NHCV

K_means(Net5[[1]], gamma[4],delta[2], 6, 5)
sort(Net5[[6]]) # Nodes with 2 overlapped membership

######################### Case 6 #########################
gamma = seq(from=0.0018,to=0.003,by=0.0003);
delta = seq(from=0.015,to=0.02,by=0.001);
Case6 <- Model_Sel(Net6[[1]],gamma,delta);

count <- 1;
MisCl_rate6 <- matrix(rep(0,length(gamma)*length(delta)),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate6[i,j] <- CV(Net6[[1]],gamma[i],delta[j],10)    
    print(count)
    count <- count + 1;
  }
}

CV_ind <- which(MisCl_rate6==min(MisCl_rate6[Case6[[4]]==3]),arr.ind=TRUE)

Eval_func(Net6,gamma[4],delta[1],6) # AIC
Eval_func(Net6,gamma[4],delta[1],6) # BIC
Eval_func(Net6,gamma[1],delta[2],6) # NHCV

K_means(Net6[[1]], gamma[1],delta[2], 5, 6)
sort(Net6[[6]]) # Nodes with 2 overlapped membership
sort(Net6[[7]]) # Nodes with 3 overlapped membership

