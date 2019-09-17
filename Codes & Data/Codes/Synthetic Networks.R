library(latex2exp)
library("plot3D")
source('ADMM_Optim.R') 

Network1 <- function(N,K,NNZ, case){
  ### Ingredients for creating network : alpha, F, D, S                      
  ### Generate alpha from uniform distribution
  alpha = runif(1,-11,-10)
  ### Create diagonal matrix D with non-negative weights
  D = diag(runif(1,19,20),K)
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
  F = (diag(N)-1/N*matrix(1,N,N))%*%F
  
  ad_hoc_ind = c();  A = choose(K,2);
  for(i in 1:(K-1)){
    for(j in i:(K-1)){
      ind = cbind( sample( ((N/K)*(i-1)+1):((N/K)*i), floor(NNZ/A), replace = FALSE),
                   sample( ((N/K)*(j)+1):((N/K)*(j+1)), floor(NNZ/A), replace = FALSE ) )
      ad_hoc_ind = rbind(ad_hoc_ind,ind)
    }
  }
  
  S = matrix(0,N,N)
  S[ad_hoc_ind]=runif(1,19,20)
  S = S + t(S)
  diag(S) = 0
  
  ### call function to generate the adjacency matrix
  X1 <- SynData(alpha,F,D,S,N)
  
  ### generate the figure to illustrate the network
  X_draw <- graph_from_adjacency_matrix(X1, mode = c("undirected"))
  
  Edge = which(S>0, arr.ind=TRUE)
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
  
  E(X_draw)$color <- "red"
  E(X_draw)[list]$color <- "green"
  
  plot(X_draw, family="serif", vertex.size=10, vertex.label=NA)
  title(sprintf("Case %d",case),family="LM Roman 10")

  Res <- list(X1,B)
  return(Res)
}

Network2<-function(N,N_K,N_K_1,K,NNZ, case){
  ### Ingredients for creating network : alpha, F, D, S                      
  ### Generate alpha from random uniform distribution
  alpha = runif(1,-11,-10)
  
  ### Create random symmetric sparse matrix S, 
  ### whose entries are greater than or equal to zero.
  ad_hoc_ind = c();  A = choose(K,2);
  for(i in 1:(K-1)){
    for(j in i:(K-1)){
      ind = cbind( sample( ((N/K)*(i-1)+1):((N/K)*i), floor(NNZ/A), replace = FALSE),
                   sample( ((N/K)*(j)+1):((N/K)*(j+1)), floor(NNZ/A), replace = FALSE ) )
      ad_hoc_ind = rbind(ad_hoc_ind,ind)
    }
  }
  
  S = matrix(0,N,N)
  S[ad_hoc_ind]=runif(1,19,20)
  S = S + t(S)
  diag(S) = 0
  
  ### Create diagonal matrix D with non-negative weights
  D = diag(runif(1,19,20),K)
  
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
  
  if(N_K_1 != 0){
    TT = sample(setdiff(1:N,ad_hoc_ind),N_K_1,replace=FALSE)
    F[TT,] = 0
    for(i in 1:length(TT)){
      F[TT[i],][c(sample(1:K,2,replace = FALSE))] = 1
    }
  }else{
    TT = 0
  }
  
  RR = sample(setdiff(1:N,c(TT,ad_hoc_ind)), N_K, replace = FALSE) 
  F[RR,] = 1
  
  F = (diag(N)-(1/N)*matrix(1,N,N))%*%F

  # S_Ind = cbind(sample(setdiff(1:N,RR),floor(length(RR)/K),replace=FALSE),
  #              sample(RR,floor(length(RR)/K),replace=FALSE))
  # S[S_Ind] = runif(1,14,15)
  # S = S + t(S)
  # diag(S) = 0
  
  ### call function to generate the adjacency matrix
  X2 <- SynData(alpha,F,D,S,N)
  ### generate the figure to illustrate the network
  X_draw2 <- graph_from_adjacency_matrix(X2, mode = c("undirected"))

  Edge = which(S>0, arr.ind=TRUE)
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
  E(X_draw2)$color <- "red"
  E(X_draw2)[list]$color <- "green"
  
  plot(X_draw2, family="serif", vertex.size=10, vertex.label=NA)
  title(sprintf("Case %d",case),family="serif")
  
  output <- list(X2,TT,RR)
  return(output)
}

Network3<-function(N,N_K,N_K_1,K,NNZ,case){
  ### Ingredients for creating network : alpha, F, D, S                      
  ### Generate alpha from random uniform distribution
  alpha = runif(1,-4,-3)
  
  ### Create random symmetric sparse matrix S, 
  ### whose entries are greater than or equal to zero.
  ad_hoc_ind = c();
  for(i in 1:(K-1)){
    for(j in i:(K-1)){
      ind = cbind( sample( ((N/K)*(i-1)+1):((N/K)*i), floor(NNZ/K), replace = FALSE),
                   sample( ((N/K)*(j)+1):((N/K)*(j+1)), floor(NNZ/K), replace = FALSE ) )
      ad_hoc_ind = rbind(ad_hoc_ind,ind)
    }
  }
  
  S = matrix(0,N,N)
  S[ad_hoc_ind]=runif(1,13,14)
  
  ### Create diagonal matrix D with non-negative weights
  D = diag(runif(1,13,14),K)
  
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
  
  if(N_K_1 != 0){
    TT = sample(setdiff(1:N,ad_hoc_ind),N_K_1,replace=FALSE)
    F[TT,] = 0
    for(i in 1:length(TT)){
      F[TT[i],][c(sample(1:K,2,replace = FALSE))] = 1
    }
  }else{
    TT = 0
  }
  
  RR = sample(setdiff(1:N,c(ad_hoc_ind,TT)), N_K, replace = FALSE) 
  F[RR,] = 1
  
  F = (diag(N)-(1/N)*matrix(1,N,N))%*%F
  
  #S_Ind = cbind(sample(setdiff(1:N,RR),floor(length(RR)/K),replace=FALSE),
  #              sample(RR,floor(length(RR)/K),replace=FALSE))
  #S[S_Ind] = runif(1,14,15)
  #S = S + t(S)
  #diag(S) = 0
  
  ### call function to generate the adjacency matrix
  X2 <- SynData(alpha,F,D,S,N)
  ### generate the figure to illustrate the network
  X_draw2 <- graph_from_adjacency_matrix(X2, mode = c("undirected"))
  
  Edge = which(S>0, arr.ind=TRUE)
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
  E(X_draw2)$color <- "red"
  E(X_draw2)[list]$color <- "green"
  #plot(X_draw2,main=TeX(sprintf("$\\textbf{Case %d : n = %d, n_{1} = %d, K = %d,\\left|S^{*}\\right| = %d}$", case, N, N, K, NNZ)))
  plot(X_draw2,main=TeX(sprintf("$\\textbf{Case  %d}$", case)))
  
  output <- list(X2,TT,RR)
  return(output)
}

SynData <- function(alpha,F,D,S,N){
  # Create empty matrix for storing adjacency matrix X and Probability for 
  X <- matrix(0,N,N); 
  P <- matrix(0,N,N);
  L <- F%*%D%*%t(F); 
  P <- exp(alpha+L+S)/(1+exp(alpha+L+S));
  
  # Create random graph according to (a, L, S)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      if(P[i,j]>runif(1,0,1))
        X[i,j]<-1
    }
  }
  X <- X + t(X)
  return(X)
}

Model_Sel <- function(X1,gamma,delta){
  lambda = 1; Count = 1; N = dim(X1)[1];
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
  
  result <- list(AIC,BIC,non_zeroS,L_Rank)
  return(result)
}

K_means<-function(X, gamma, delta, K_clusters, case){
  result <- ADMM(X, gamma, delta)
  a<-result[[1]]
  M<-result[[2]]
  L<-result[[3]]
  S<-result[[4]]
  
  K = qr(L)$rank
  N = ncol(X)
  FAC = eigen(L)$vectors[,1:K]
  pca <- prcomp(FAC,scale=TRUE)

  if(case==1 | case == 2 | case ==3){
    KMeans = kmeans(FAC, K_clusters, iter.max = 1000, nstart = 100, 
                    algorithm = "Hartigan-Wong")
    
    for(i in 1:K_clusters){
      print(length(which(KMeans$cluster==i,arr.ind=TRUE)))
      print(which(KMeans$cluster==i,arr.ind=TRUE))
    }
    
    plot(FAC[,1], FAC[,2], family="serif", col=KMeans$cluster, cex=1, pch=0, xlab="PC1", ylab="PC2")
    title(sprintf("Case %d",case),family="serif")
  }
  else{
    KMeans = kmeans(cbind(pca$x[,1],pca$x[,2]), K_clusters, iter.max = 1000, nstart = 100, 
                    algorithm = "Hartigan-Wong")
    
    for(i in 1:K_clusters){
      print(length(which(KMeans$cluster==i,arr.ind=TRUE)))
      print(which(KMeans$cluster==i,arr.ind=TRUE))
    }
    
    plot(pca$x[,1], pca$x[,2],family="serif", col=KMeans$cluster, cex=1, pch=0, xlab="PC1", ylab="PC2")
    title(sprintf("Case %d",case),family="serif")
  }
  #scatter3D(eigen(L)$vectors[,1], eigen(L)$vectors[,2], eigen(L)$vectors[,3], theta = 15, phi = 20)
}


Eval <- function(X,Result,gamma,delta,K_Scree,low,upp){
  gamma = seq(from=0.01,to=0.02,by=0.001);
  delta = seq(from=0.02,to=0.03,by=0.001);
  
  X <- Net9[[1]]; Ad_hoc_true <- Net9[[4]]; Result <- Res9;
  AIC_indx <- which(Result[[1]] == min(Result[[1]]),arr.ind=TRUE);
  BIC_indx <- which(Result[[2]] == min(Result[[2]]),arr.ind=TRUE);
  Heu_indx <- which(Result[[4]] == K_Scree,arr.ind=TRUE)
  a = (Result[[3]]>=low*(sum(X)/2)); b = (Result[[3]]<=upp*(sum(X)/2)); c = (Result[[4]]==K_Scree);
  d = (a == b); which(d==c,arr.ind=TRUE);
  
  ### AIC ###
  Res <- ADMM(X,gamma[AIC_indx[1]],delta[AIC_indx[2]]);
  L = Res[[3]]; S = Res[[4]];
  AIC_L_Rank <- qr(L)$rank;
  est_S_AIC = which(S!=0,arr.ind=TRUE)
  
  ### BIC ###
  gamma = BIC_indx[1]; delta = BIC_indx[2];  
  Res <- ADMM(X,gamma,delta);
  L = Res[[3]]; S = Res[[4]];
  BIC_L_Rank <- qr(L)$rank;
  est_S_BIC = which(S!=0,arr.ind=TRUE)
  
  ### Heuristic ###
}

Eval <- function(Net,Res,gamma,delta){
  Net = Net3; Res = Res3;
  X <- Net[[1]]; Ad_hoc_true <- Net[[2]]; Result <- Res; N = ncol(X);
  AIC_indx <- which(Result[[1]] == min(Result[[1]]),arr.ind=TRUE);
  BIC_indx <- which(Result[[2]] == min(Result[[2]]),arr.ind=TRUE);
  
  ### AIC ###
  Res <- ADMM(X,gamma[AIC_indx[1]],delta[AIC_indx[2]]);
  L = Res[[3]]; S = Res[[4]];
  AIC_L_Rank <- qr(L)$rank;
  est_S_AIC = which(S!=0,arr.ind=TRUE)
  Ad_hoc_est <- est_S_AIC[which(est_S_AIC[,1]<est_S_AIC[,2],arr.ind=TRUE),]
  M2_AIC <- ((length(Ad_hoc_true) - length(setdiff(Ad_hoc_true,Ad_hoc_est)))/2)/(length(Ad_hoc_true)/2)
  #M3_AIC <- length()/((N*(N-1)/2)-(length(Ad_hoc_true)/2))
  AIC <- c(AIC_L_Rank,M2_AIC);
  
  ### BIC ###
  Res <- ADMM(X,gamma[BIC_indx[1]],delta[BIC_indx[2]]);
  L = Res[[3]]; S = Res[[4]];
  BIC_L_Rank <- qr(L)$rank;
  est_S_BIC = which(S!=0,arr.ind=TRUE)
  Ad_hoc_est <- est_S_BIC[which(est_S_BIC[,1]<est_S_BIC[,2],arr.ind=TRUE),]
  M2_BIC <- ((length(Ad_hoc_true) - length(setdiff(Ad_hoc_true,Ad_hoc_est)))/2)/(length(Ad_hoc_true)/2)
  #M3_BIC <- length()/((N*(N-1)/2)-(length(Ad_hoc_true)/2))
  
  BIC <- c(BIC_L_Rank,M2_BIC);
  RESULT<-list(AIC,BIC)
  return(RESULT)
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

Intersected_set <- function(d,c){
  f = matrix(0,nrow = length(gamma), ncol = length(delta))
  for(i in 1:length(gamma)){
    for(j in 1:length(delta)){
      if(d[i,j]==1 & c[i,j]==1){
        f[i,j] = 1
      }else{
        f[i,j] = 0
      }
    }
  }
  return(f)
}

















