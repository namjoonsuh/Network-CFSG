# Citation Network
## Introduction 
This is for short descriptions on running codes for "Numerical analysis and application" section in our manuscript paper.

## Synthetic Data
We describe a code for running an experiment on synthetic data, where n=80, K=4, nnz=20 (number of non-zero elements in upper-triangular part of S). Our goal in this section is to check our model can capture two disctint dependent structure in network, and correctly classify each nodes to correct clusters.
First, we need to load following library and functions for running codes
```R
library('matrix')       ## For creating random sparse matrix S
library('igraph')       ## For converting adjacency matrix into network object and for drawing graph
source('SynData.R')     ## Function for creating binary adjacency matrix
source('ADMM_Optim.R')  ## Main function of ADMM algorithm for estimation
source('functions.R')   ## Function for drawing screeplot 
```

```R
#################### data generation N = 80, K = 4, NNZ = 20 ####################
N = 80; K = 4; NNZ = 20;
### Ingredients for creating network : alpha, F, D, S                      
### Generate alpha from random uniform distribution
alpha = runif(1,-3,-2)
### Create diagonal matrix D with non-negative weights
D = diag(runif(1,7,8),K)
### Create binary matrix F(0,1) whose columns are orthogonal with each other.
C = floor(N/K);               # C : number of 1s for each column
F = matrix(0,N,K);            # Create empty F matrix
if(K>0){
  for(l in 1:K){
    F[((l-1)*C+1):(l*C),l] <- 1
  }
  if(N%%C!=0)
    F[(K*C+1):N,K] <- 1       # In last column of F, we fill the last remaining entries with 1s, 
  # when N is not divided by K.
}
if(K==0){
  F <- matrix(0,N,K)
}
RR = sample(1:N, 1, replace = FALSE) 
F[RR,] = 1
F = (diag(N)-1/N*matrix(1,N,N))%*%F
### Create random symmetric sparse matrix S ### 
S = matrix(0,N,N)
S = rsparsematrix(N, N, nnz = NNZ, symmetric = TRUE, rand.x = runif)
S[which(S!=0,arr.ind=TRUE)] = runif(1,6,7)
diag(S) = 0
### call function to generate the adjacency matrix
X2 <- SynData(alpha,F,D,S,N)
### generate the figure to illustrate the network
X_draw2 <- graph_from_adjacency_matrix(X2, mode = c("undirected"))
plot(X_draw2)
```
Then we draw the screeplot to get the knack of how many topics there are, and choose a proper pair of tuning parameter.
This case (gamma,delta) = (0.004, 0.014). Here we omit the code for choosing tuning parameter. (Detailed code is in https://github.com/namjoonsuh/Citation-Network/blob/master/Codes%20%26%20Data/Codes/Synthetic%20Data.R)

```R
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

#### Perform K-means clustering for assigning the membership
KMeans = kmeans(eigen(L2)$vectors[,1:K], K, iter.max = 1000, nstart = 100, algorithm = "Hartigan-Wong")
plot(eigen(L2)$vectors[,1:K], col=KMeans$cluster,pch=2, main="N=80,K=4,NNZ=20",
     xlab="First_eigen vector", ylab="Second_eigen vector")

################## Color the estimated ad-hoc edges with blue ################## 
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
```
Following three plots are in sequence, screeplot of eigenvalues on adjacency matrix X2, a plot of result on K-means clustering on first K(=rank(estimated L)) eigenvectors of estimated L matrix. Lastly, We color the blue on edges whose corresponding entries of estimated S matrix are positive.

![](https://github.com/namjoonsuh/Citation-Network/blob/master/Manuscript%20Tex%20file/Figures/Fig_git_1.png)
![](https://github.com/namjoonsuh/Citation-Network/blob/master/Manuscript%20Tex%20file/Figures/Fig_git_3.png)
![](https://github.com/namjoonsuh/Citation-Network/blob/master/Manuscript%20Tex%20file/Figures/Fig_git_2.png)


## Data reference for section 6.2.
Ji, P. and Jin, J. (2016). Coauthorship and citation networks for statisticians. Annals of Applied Statistics with discussion, Vol. 10, No. 4, 1779-1812. Data set available from :
http://faculty.franklin.uga.edu/psji/statistician-coauthorship-and-citation
