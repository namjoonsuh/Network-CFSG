####################  preparation  ####################
rm(list=ls())
setwd('C:/Users/nsuh3/Desktop/Citation Network/FLag/Codes & Data/Codes')
library('Matrix')      ## For matrix computation on BIC
library('igraph')      ## For drawing graph
source('SynData.R')    ## Function for creating binary network
source('ADMM_Optim.R') ## Function of ADMM algorithm for estimation
source('functions.R')
#################### data generation N = 30, K = 3, NNZ = 10 ####################
set.seed(1234)
par(mfrow=c(1,3))
N = 30; K = 3; NNZ = 10;
### Ingredients for creating network : alpha, F, D, S                      
### Generate alpha from random uniform distribution
alpha = runif(1,-3,-2)
### Create diagonal matrix D with non-negative weights
D = diag(runif(1,7,8),K)
### Create binary matrix F(0,1) whose columns are orthogonal with each other.
C = floor(N/K);              # C : number of 1s for each column
F = matrix(0,N,K);           # Create empty F matrix
if(K>0){
  for(l in 1:K){
    F[((l-1)*C+1):(l*C),l] <- 1
  }
  if(N%%C!=0)
    F[(K*C+1):N,K] <- 1      # In last column of F, we fill the last remaining entries with 1s, 
  # when N is not divided by K.
}
if(K==0){
  F <- matrix(0,N,K)
}
RR = sample(1:N, 1, replace = FALSE) 
F[RR,] = 1
F = (diag(N)-1/N*matrix(1,N,N))%*%F
### Create random symmetric sparse matrix S, whose entries are greater than or equal to zero.
S = matrix(0,N,N)
S = rsparsematrix(N, N, nnz = NNZ, symmetric = TRUE, rand.x = runif)
S[which(S!=0,arr.ind=TRUE)] = runif(1,6,7)
diag(S) = 0
### call function to generate the adjacency matrix
X1 <- SynData(alpha,F,D,S,N)
### generate the figure to illustrate the network
X_draw <- graph_from_adjacency_matrix(X1, mode = c("undirected"))
plot(X_draw)

#################### Choose proper tuning parameter #########################
plotScree(X1)
gamma = seq(from=0.01,to=0.1,by=0.01);
delta = seq(from=0.01,to=0.1,by=0.01);
lambda = 1; Count = 1;
AIC <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
BIC <- matrix( 0, nrow=length(gamma),ncol=length(delta) )

non_zeroS <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
L_Rank <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
Like <- matrix( 0, nrow=length(gamma),ncol=length(delta) )

for(g in 1:length(gamma)){
  for(d in 1:length(delta)) {
    ### Use the ADMM method to estimate the parameters ###
    result <- ADMM(X1, gamma[g], delta[d])
    
    a<-result[[1]]
    M<-result[[2]]
    L<-result[[3]]
    S<-result[[4]]
    
    non_zero_S <- 0
    for(i in 1:N){
      for(j in 1:N){
        if(i<j){
          if(S[i,j] != 0)
            non_zero_S = non_zero_S + 1
        }
      }
    }
    
    K<-qr(L)$rank
    M_absolute <- (N*K-((K-1)*K/2)) + non_zero_S + 1
    
    log_sum<-log(1+exp(a + L + S))
    
    log_max <- 0;
    log_max <- a*sum(X1[upper.tri(X1,diag=FALSE)]) +
      (1/2)*sum(X1*(L+S)) - sum(log_sum[upper.tri(log_sum,diag=FALSE)])
    
    AIC[g,d] <- -2*log_max + M_absolute*2
    BIC[g,d] <- -2*log_max + M_absolute*(log((N*(N-1))/2))
    Like[g,d] <- log_max
    non_zeroS[g,d] <- non_zero_S
    L_Rank[g,d] <- qr(L)$rank
    
    print(Count)
    Count <- Count + 1
  }
}
################## Assign memberships to each nodes ##################
result1 <- ADMM(X1, 0.01, 0.03) ## gamma : 0.01, delta : 0.03 ##
a1<-result1[[1]]
M1<-result1[[2]]
L1<-result1[[3]]
S1<-result1[[4]]

K = qr(L1)$rank

N = ncol(X1)
First_eigen1 = eigen(L1)$vectors[,1]
Second_eigen1 = eigen(L1)$vectors[,2]
Third_eigen1 = eigen(L1)$vectors[,3]
Label1 = seq(from=1,to=N,by=1)

KMeans = kmeans(eigen(L1)$vectors[,1:K], K, iter.max = 1000, nstart = 100, algorithm = "Hartigan-Wong")
plot(eigen(L1)$vectors[,1:K], col=KMeans$cluster,pch=2, main="N=30,K=3,NNZ=10",
     xlab="First_eigen vector", ylab="Second_eigen vector")

################## Draw the plot with ad-hoc edges ################## 
Edge = which(S1>0, arr.ind=TRUE)
g = graph_from_edgelist(Edge)
S_g = as.undirected(g)

A=get.edgelist(X_draw)
B=get.edgelist(S_g)

list = rep(0,length(B[,1]))

for(i in 1:length(B[,1])){
  for(j in 1:length(A[,1])){
    if(B[i,1]==A[j,1] && B[i,2]==A[j,2]){
      list[i] = j;
      break;
    }
  }
}
X_draw <- graph_from_adjacency_matrix(X1, mode = c("undirected"))
E(X_draw)$color <- "black"
E(X_draw)[list]$color <- "blue"
par(mfrow=c(1,3))
plot(X_draw,main="N=30,K=3,NNZ=10")




#################### data generation N = 80, K = 4, NNZ = 20 ####################
N = 80; K = 4; NNZ = 20;
### Ingredients for creating network : alpha, F, D, S                      
### Generate alpha from random uniform distribution
alpha = runif(1,-3,-2)
### Create diagonal matrix D with non-negative weights
D = diag(runif(1,7,8),K)
### Create binary matrix F(0,1) whose columns are orthogonal with each other.
C = floor(N/K);              # C : number of 1s for each column
F = matrix(0,N,K);           # Create empty F matrix
if(K>0){
  for(l in 1:K){
    F[((l-1)*C+1):(l*C),l] <- 1
  }
  if(N%%C!=0)
    F[(K*C+1):N,K] <- 1      # In last column of F, we fill the last remaining entries with 1s, 
  # when N is not divided by K.
}
if(K==0){
  F <- matrix(0,N,K)
}
RR = sample(1:N, 1, replace = FALSE) 
F[RR,] = 1
F = (diag(N)-1/N*matrix(1,N,N))%*%F
### Create random symmetric sparse matrix S, whose entries are greater than or equal to zero.
S = matrix(0,N,N)
S = rsparsematrix(N, N, nnz = NNZ, symmetric = TRUE, rand.x = runif)
S[which(S!=0,arr.ind=TRUE)] = runif(1,6,7)
diag(S) = 0
### call function to generate the adjacency matrix
X2 <- SynData(alpha,F,D,S,N)
### generate the figure to illustrate the network
X_draw2 <- graph_from_adjacency_matrix(X2, mode = c("undirected"))
plot(X_draw2)

#################### Choose proper tuning parameter #########################
plotScree(X2)
gamma = seq(from=0.001,to=0.01,by=0.001);
delta = seq(from=0.01,to=0.02,by=0.001);
lambda = 1; Count = 1;
AIC2 <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
BIC2 <- matrix( 0, nrow=length(gamma),ncol=length(delta) )

non_zeroS2 <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
L_Rank2 <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
Like2 <- matrix( 0, nrow=length(gamma),ncol=length(delta) )

for(g in 1:length(gamma)){
  for(d in 1:length(delta)) {
    ### Use the ADMM method to estimate the parameters ###
    result <- ADMM(X2, gamma[g], delta[d])
    
    a<-result[[1]]
    M<-result[[2]]
    L<-result[[3]]
    S<-result[[4]]
    
    non_zero_S <- 0
    for(i in 1:N){
      for(j in 1:N){
        if(i<j){
          if(S[i,j] != 0)
            non_zero_S = non_zero_S + 1
        }
      }
    }
    
    K<-qr(L)$rank
    M_absolute <- (N*K-((K-1)*K/2)) + non_zero_S + 1
    
    log_sum<-log(1+exp(a + L + S))
    
    log_max <- 0;
    log_max <- a*sum(X2[upper.tri(X1,diag=FALSE)]) +
      (1/2)*sum(X2*(L+S)) - sum(log_sum[upper.tri(log_sum,diag=FALSE)])
    
    AIC2[g,d] <- -2*log_max + M_absolute*2
    BIC2[g,d] <- -2*log_max + M_absolute*(log((N*(N-1))/2))
    Like2[g,d] <- log_max
    non_zeroS2[g,d] <- non_zero_S
    L_Rank2[g,d] <- qr(L)$rank
    
    print(Count)
    Count <- Count + 1
  }
}

################## Assign memberships to each nodes ##################
plotScree(X2)
result2 <- ADMM(X2, 0.004, 0.014) ## gamma : 0.004, delta : 0.014 ##
a2<-result2[[1]]
M2<-result2[[2]]
L2<-result2[[3]]
S2<-result2[[4]]

K = qr(L2)$rank

N = ncol(X2)
First_eigen1 = eigen(L2)$vectors[,1]
Second_eigen1 = eigen(L2)$vectors[,2]
Third_eigen1 = eigen(L2)$vectors[,3]
Label1 = seq(from=1,to=N,by=1)

KMeans = kmeans(eigen(L2)$vectors[,1:K], K, iter.max = 1000, nstart = 100, algorithm = "Hartigan-Wong")
plot(eigen(L2)$vectors[,1:K], col=KMeans$cluster,pch=2, main="N=80,K=4,NNZ=20",
     xlab="First_eigen vector", ylab="Second_eigen vector")

################## Draw the plot with ad-hoc edges ################## 
Edge = which(S2>0, arr.ind=TRUE)
g = graph_from_edgelist(Edge)
S_g = as.undirected(g)

A=get.edgelist(X_draw2)
B=get.edgelist(S_g)

list = rep(0,length(B[,1]))

for(i in 1:length(B[,1])){
  for(j in 1:length(A[,1])){
    if(B[i,1]==A[j,1] && B[i,2]==A[j,2]){
      list[i] = j;
      break;
    }
  }
}
E(X_draw2)$color <- "black"
E(X_draw2)[list]$color <- "blue"
plot(X_draw2,main="N=80,K=4,NNZ=20")





#################### data generation N = 120, K = 5, NNZ = 40 ###########################
set.seed(1234)
N = 120; K = 5; NNZ = 40;
### Ingredients for creating network : alpha, F, D, S                      
### Generate alpha from random uniform distribution
alpha = runif(1,-3,-2)
### Create diagonal matrix D with non-negative weights
D = diag(runif(1,7,8),K)
### Create binary matrix F(0,1) whose columns are orthogonal with each other.
C = floor(N/K);              # C : number of 1s for each column
F = matrix(0,N,K);           # Create empty F matrix
if(K>0){
  for(l in 1:K){
    F[((l-1)*C+1):(l*C),l] <- 1
  }
  if(N%%C!=0)
    F[(K*C+1):N,K] <- 1      # In last column of F, we fill the last remaining entries with 1s, 
  # when N is not divided by K.
}
if(K==0){
  F <- matrix(0,N,K)
}
RR = sample(1:N, 1, replace = FALSE) 
F[RR,] = 1
F = (diag(N)-1/N*matrix(1,N,N))%*%F
### Create random symmetric sparse matrix S, whose entries are greater than or equal to zero.
S = matrix(0,N,N)
S = rsparsematrix(N, N, nnz = NNZ, symmetric = TRUE, rand.x = runif)
S[which(S!=0,arr.ind=TRUE)] = runif(1,6,7)
diag(S) = 0
### call function to generate the adjacency matrix
X3 <- SynData(alpha,F,D,S,N)
### generate the figure to illustrate the network
X_draw3 <- graph_from_adjacency_matrix(X3, mode = c("undirected"))
plot(X_draw3)

#################### Choose proper tuning parameter #########################
plotScree(X3)
gamma = seq(from=0.001,to=0.01,by=0.001);
delta = seq(from=0.01,to=0.02,by=0.001);
lambda = 1; Count = 1;
AIC3 <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
BIC3 <- matrix( 0, nrow=length(gamma),ncol=length(delta) )

non_zeroS3 <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
L_Rank3 <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
Like3 <- matrix( 0, nrow=length(gamma),ncol=length(delta) )

for(g in 1:length(gamma)){
  for(d in 1:length(delta)) {
    ### Use the ADMM method to estimate the parameters ###
    result <- ADMM(X3, gamma[g], delta[d])
    
    a<-result[[1]]
    M<-result[[2]]
    L<-result[[3]]
    S<-result[[4]]
    
    non_zero_S <- 0
    for(i in 1:N){
      for(j in 1:N){
        if(i<j){
          if(S[i,j] != 0)
            non_zero_S = non_zero_S + 1
        }
      }
    }
    
    K<-qr(L)$rank
    M_absolute <- (N*K-((K-1)*K/2)) + non_zero_S + 1
    
    log_sum<-log(1+exp(a + L + S))
    
    log_max <- 0;
    log_max <- a*sum(X3[upper.tri(X1,diag=FALSE)]) +
      (1/2)*sum(X3*(L+S)) - sum(log_sum[upper.tri(log_sum,diag=FALSE)])
    
    AIC3[g,d] <- -2*log_max + M_absolute*2
    BIC3[g,d] <- -2*log_max + M_absolute*(log((N*(N-1))/2))
    Like3[g,d] <- log_max
    non_zeroS3[g,d] <- non_zero_S
    L_Rank3[g,d] <- qr(L)$rank
    
    print(Count)
    Count <- Count + 1
  }
}
#################### Assign membership to each nodes ####################
result3 <- ADMM(X3, 0.003, 0.013) ## gamma : 0.003, delta : 0.013 ##
a3<-result3[[1]]
M3<-result3[[2]]
L3<-result3[[3]]
S3<-result3[[4]]

K = qr(L3)$rank

N = ncol(X3)
First_eigen1 = eigen(L3)$vectors[,1]
Second_eigen1 = eigen(L3)$vectors[,2]
Third_eigen1 = eigen(L3)$vectors[,3]
Label3 = seq(from=1,to=N,by=1)

KMeans = kmeans(eigen(L3)$vectors[,1:K], K, iter.max = 1000, nstart = 100, algorithm = "Hartigan-Wong")
plot(eigen(L3)$vectors[,1:K], col=KMeans$cluster,pch=2, main="N=120, K=5, NNZ=40",
     xlab="First_eigen vector", ylab="Second_eigen vector")

################## Draw the plot with ad-hoc edges ################## 
Edge = which(S3>0, arr.ind=TRUE)
g = graph_from_edgelist(Edge)
S_g = as.undirected(g)

#Ad_hoc_S <- graph_from_edgelist(el, directed = TRUE)
#plot(Ad_hoc_S)

A=get.edgelist(X_draw3)
B=get.edgelist(S_g)

list = rep(0,length(B[,1]))

for(i in 1:length(B[,1])){
  for(j in 1:length(A[,1])){
    if(B[i,1]==A[j,1] && B[i,2]==A[j,2]){
      list[i] = j;
      break;
    }
  }
}
X_draw3 <- graph_from_adjacency_matrix(X3, mode = c("undirected"))
E(X_draw3)$color <- "black"
E(X_draw3)[list]$color <- "blue"
plot(X_draw3, main = "N=120, K=5, NNZ=120")

par(mfrow=c(1,1))
