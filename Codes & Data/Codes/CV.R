source('functions.R')
CV <- function(X,gamma,delt,K){
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
gamma = seq(from=0.01,to=0.02,by=0.001);
delta = seq(from=0.02,to=0.03,by=0.001);

count <- 1;
MisCl_rate <- matrix(rep(0,length(gamma)^2),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate[i,j] <- CV(Net1[[1]],gamma[i],delta[j],20)    
    print(count)
    count <- count + 1;
  }
}
Min_case1 <- which(MisCl_rate==min(MisCl_rate),arr.ind=TRUE) # gamma : 0.011, delta : 0.021
C1 <- ADMM(Net1[[1]],0.011,0.021)
qr(C1[[3]])$rank; C1[[4]] 

######################### Scenario 2 #########################
gamma = seq(from=0.004,to=0.005,by=0.0001)
delta = seq(from=0.009,to=0.01,by=0.0001)

count <- 1;
MisCl_rate <- matrix(rep(0,length(gamma)^2),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate[i,j] <- CV(Net2[[1]],gamma[i],delta[j],20)    
    print(count)
    count <- count + 1;
  }
}


