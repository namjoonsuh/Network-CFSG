rm(list=ls())
setwd('C:/Users/Namjoon Suh/Desktop/Stat+PDE/Citation-Network/Codes & Data/Codes')
source('ADMM_Optim.R')
source('functions.R')
library("igraph"); library("CINNA");

pol = as.matrix(read.table("polblogsgiant.txt"));
gamma = seq(from=0.00061,to=0.0007,by=0.00001);
delta = seq(from=0.00651,to=0.0075,by=0.0001);
Res_pol <- Model_Sel(pol,gamma,delta);
save(Res_pol, file = "PolBlog_result.RData")
load("PolBlog_result.RData") 


Mat = as_adjacency_matrix(gian_pol[[1]],type="both")
M = as.matrix(Mat)

for(i in 1:1222){
  for(j in 1:1222){
    if(M[i,j]==1){
      M[j,i] = 1;
    }
  }
} 

Res = ADMM(M,0.00065,0.005)
qr(Res[[3]])$rank
count <- 0;
for(i in 1:1222){
  for(j in 1:1222){
    if(i<j){
      if(Res[[4]][i,j]!=0){
        count <- count + 1
      }
    }  
  }
}
count # Number of Ad-hoc edges 160 

PC1 <- eigen(Res[[3]])$vectors[,1];
PC2 <- eigen(Res[[3]])$vectors[,2];
A = cbind(PC1,PC2);
plot(PC1,PC2);

# Pairs of blogs 
Edge = which(Res[[4]]>0, arr.ind=TRUE)
g = graph_from_edgelist(Edge)
S_g = as.undirected(g)
Edge_list = get.edgelist(S_g)

Ad_hoc_link <- data.frame(Pair1 = Edge_list[,1], Pair2 = Edge_list[,2],  S_val = Res[[4]][Edge_list])
Sorted_Ad_hoc_link <- Ad_hoc_link[order(-Ad_hoc_link$S_val),] 
Edge_list <- as.matrix(Sorted_Ad_hoc_link[,1:2])
table(Edge_list)


KMeans = kmeans(A, 2, iter.max = 1000, nstart = 100, algorithm = "Hartigan-Wong")
Topic1 = which(KMeans$cluster == 1,arr.ind = TRUE)
Topic2 = which(KMeans$cluster == 2,arr.ind = TRUE)



