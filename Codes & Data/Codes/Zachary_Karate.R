rm(list=ls())
library(igraphdata) # for karate data 
library(igraph)
source('Synthetic Networks.R')
source('ADMM_Optim.R') # delta = 10^-7, r = 15
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

##### Draw graph ######
##### adj : adjacency matrix, S : Sparse ad-hoc matrix
##### C1 : Nodes in cluster 1, C2 : Nodes in cluster 2
plot_kara <- function(adj,S,C1,C2){
  rownames(adj)<-seq(1:34)
  colnames(adj)<-seq(1:34)
  
  Kar_graph <- graph_from_adjacency_matrix(adj, mode = c("undirected"))
  
  Edge = which(S>0, arr.ind=TRUE)
  g = graph_from_edgelist(Edge)
  S_g = as.undirected(g)
  
  A=get.edgelist(Kar_graph)
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
  
  V(Kar_graph)[C1]$color <- "blue"
  V(Kar_graph)[C2]$color <- "orange"
  
  E(Kar_graph)$color <- "gray"
  E(Kar_graph)[list]$color <- "green"
  
  plot(Kar_graph,edge.width=2)
}

##### Load the Karate dataset #####
data(karate) # type: ?karate to see description 
Karate_ad <- as.matrix(as_adjacency_matrix(karate, type = c("both", "upper", "lower"),
                    attr = NULL, edges = FALSE, names = TRUE,
                    sparse = igraph_opt("sparsematrices")))

plotScree(Karate_ad, 1)

##### Model Selection through HNCV #####
gamma = seq(from=0.012,to=0.0128,by=0.0002);
delta = seq(from=0.04,to=0.05,by=0.002);

Kar_model <- Model_Sel(Karate_ad, gamma, delta);

count <- 1;
MisCl_rate1 <- matrix(rep(0,length(gamma)*length(delta)),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 5:length(delta)){
    MisCl_rate1[i,j] <- CV(Karate_ad,gamma[i],delta[j],10)
    print(count)
    count <- count + 1;
  }
}

CV_ind <- which(MisCl_rate1==min(MisCl_rate1[Kar_model[[4]][1:5,]==2]),arr.ind=TRUE)
Kara <- ADMM(Karate_ad,0.0126,0.048)

##### discovered ad-hoc list
list <- c();
for(i in 1:34){
  for(j in 1:34){
    if(i<j && Kara[[4]][i,j]!=0){
      list = rbind(list, c(i,j));
    }
  }
}
print(list)

##### Cluster nodes through K-means clustering #####
K=qr(Kara[[3]])$rank
vc <- eigen(Kara[[3]])$vectors[,1:2]
plot(vc[,1],vc[,2])

a<-kmeans(vc[,1:2], 2, iter.max = 1000, nstart = 100, 
       algorithm = "Hartigan-Wong")$cluster

C1 <- which(a==1,arr.ind=TRUE)
C2 <- which(a==2,arr.ind=TRUE)

plot(vc[,1], vc[,2],family="serif", col=a$cluster, cex=1, pch=0, xlab="PC1", ylab="PC2")

##### Plot Final graph of the selected model #####
plot_kara(Karate_ad,Kara[[4]],C1,C2)

