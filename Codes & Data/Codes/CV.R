source('ADMM_Optim.R')
source('functions.R')
source('Synthetic Networks.R')
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

######################### Scenario 1 #########################
######################### Case 1 #########################
gamma = seq(from=0.01,to=0.02,by=0.001);
delta = seq(from=0.02,to=0.03,by=0.001);

Case1 <- Model_Sel(Net1[[1]],gamma,delta);
Eval_func(Net1,gamma[Case1[[1]][1,1]],delta[Case1[[1]][1,2]],3,1) # AIC
Eval_func(Net1,gamma[Case1[[2]][1,1]],delta[Case1[[2]][1,2]],3,1) # BIC
Eval_func(Net1,gamma[4],delta[9],3,1) # Heuristic

count <- 1;
MisCl_rate1 <- matrix(rep(0,length(gamma)^2),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate1[i,j] <- CV(Net1[[1]],gamma[i],delta[j],10)    
    print(count)
    count <- count + 1;
  }
}

CV_ind <- which(MisCl_rate1==min(MisCl_rate1),arr.ind=TRUE)
Res1 <- ADMM(Net1[[1]],gamma[CV_ind[1,1]],delta[CV_ind[1,2]])
qr(Res1[[3]])$rank; 
Eval_func(Net1,gamma[CV_ind[1,1]],delta[CV_ind[1,2]],3,1) # Heuristic
 
######################### Case 2 #########################
gamma = seq(from=0.005,to=0.006,by=0.0001)
delta = seq(from=0.009,to=0.01,by=0.0001)

Case2 <- Model_Sel(Net2[[1]],gamma,delta);
Eval_func(Net2,gamma[Case2[[1]][1,1]],delta[Case2[[1]][1,2]],1) # AIC
Eval_func(Net2,gamma[Case2[[2]][1,1]],delta[Case2[[2]][1,2]],1) # BIC
Eval_func(Net1,gamma[4],delta[9],1) # Heuristic

count <- 1;
MisCl_rate2 <- matrix(rep(0,length(gamma)^2),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate2[i,j] <- CV(Net2[[1]],gamma[i],delta[j],10)    
    print(count)
    count <- count + 1;
  }
}

CV_ind <- which(MisCl_rate2==min(MisCl_rate2),arr.ind=TRUE)
Res2 <- ADMM(Net2[[1]],gamma[11],delta[11])
Eval_func(Net2,gamma[11],delta[11],1) # Heuristic

######################### Case 3 #########################
gamma = seq(from=0.003,to=0.004,by=0.0001)
delta = seq(from=0.0044,to=0.0046,by=0.00002)

Case3 <- Model_Sel(Net3[[1]],gamma,delta);

Eval_func(Net3,gamma[Case3[[1]][1,1]],delta[Case3[[1]][1,2]],1) # AIC
Eval_func(Net3,gamma[Case3[[2]][1,1]],delta[Case3[[2]][1,2]],1) # BIC
Eval_func(Net1,gamma[4],delta[9],1) # Heuristic

count <- 1;
MisCl_rate3 <- matrix(rep(0,length(gamma)^2),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate3[i,j] <- CV(Net3[[1]],gamma[i],delta[j],10)    
    print(count)
    count <- count + 1;
  }
}

CV_ind <- which(MisCl_rate3==min(MisCl_rate3),arr.ind=TRUE)
Eval_func(Net3,gamma[10],delta[5],1) # Heuristic


######################### Scenario 2 #########################
######################### Case 4 #########################
gamma = seq(from=0.002,to=0.004,by=0.0002);
delta = seq(from=0.015,to=0.02,by=0.0005);
Case4 <- Model_Sel(Net4[[1]][[1]],gamma,delta);
Case4[[3]]
Eval_func(Net4[[1]],gamma[Case4[[1]][1,1]],delta[Case4[[1]][1,2]],1) # AIC
Eval_func(Net4[[1]],gamma[Case4[[2]][1,1]],delta[Case4[[2]][1,2]],1) # BIC
Net4








count <- 1;
MisCl_rate4 <- matrix(rep(0,length(gamma)^2),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate4[i,j] <- CV(Net4[[1]][[1]],gamma[i],delta[j],10)    
    print(count)
    count <- count + 1;
  }
}

which(MisCl_rate4==min(MisCl_rate4),arr.ind=TRUE)
Res4 <- ADMM(Net4[[1]][[1]],gamma[[5]],delta[[9]])
qr(Res4[[3]])$rank; 
which(Res4[[4]]!=0,arr.ind=TRUE)

######################### Case 5 #########################
gamma = seq(from=0.00204,to=0.00244,by=0.00004); 
delta = seq(from=0.045,to=0.05,by=0.0005);

count <- 1;
MisCl_rate5 <- matrix(rep(0,length(gamma)^2),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate5[i,j] <- CV(Net5[[1]][[1]],gamma[i],delta[j],10)    
    print(count)
    count <- count + 1;
  }
}

## Network CV
which(MisCl_rate5==min(MisCl_rate5),arr.ind=TRUE)
Res5 <- ADMM(Net5[[1]][[1]],gamma[6],delta[2])
qr(Res5[[3]])$rank; 
K_means(Net5[[1]][[1]],gamma[6],delta[2],3,1); 
which(Res5[[4]]!=0,arr.ind=TRUE)

qr(Res2[[3]])$rank;
K_means(Net2[[1]],gamma[7],delta[11],3,1); 
which(Res2[[4]]!=0,arr.ind=TRUE)



## Heuristic
Res5 <- ADMM(Net5[[1]][[1]],0.00212,0.048)
qr(Res5[[3]])$rank; 
which(Res5[[4]]!=0,arr.ind=TRUE)

######################### Case 6 #########################
gamma = seq(from=0.002,to=0.004,by=0.0002);
delta = seq(from=0.015,to=0.02,by=0.0005);

count <- 1;
MisCl_rate6 <- matrix(rep(0,length(gamma)^2),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate6[i,j] <- CV(Net6[[1]][[1]],gamma[i],delta[j],20)    
    print(count)
    count <- count + 1;
  }
}





