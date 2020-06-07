# Karate club data example
- We take a simple example on the application of our model to famous [Zachary's Karate club dataset](https://en.wikipedia.org/wiki/Zachary%27s_karate_club). First we import necessary libraries and functions for the analysis, load the network data, and make the adjacency matrix from the network data. 
```R
##### Load libraries and functions for analysis 
library(igraphdata) # for karate data 
library(igraph)
source('Synthetic Networks.R')
source('ADMM_Optim.R') # delta = 10^-7, r = 15

##### Load the Karate dataset 
data(karate) # type: ?karate to see description 
Karate_ad <- as.matrix(as_adjacency_matrix(karate, type = c("both", "upper", "lower"),
                    attr = NULL, edges = FALSE, names = TRUE,
                    sparse = igraph_opt("sparsematrices")))
```

- Over the given grid points with gamma being ranged from 0.012 to 0.0128 with 0.0002 interval and delta being ranged from 0.04 to 0.05 with 0.002 interval, we perform a model selection through Heuristic Network Cross-validation. (There are 30 points in total on the grid.)
    1. Run a [Model_Sel](https://github.com/namjoonsuh/Citation-Network/blob/master/Codes%20%26%20Data/Codes/Synthetic%20Networks.R) code for recording the rank of the estimated L and the number of non-zero entries of the estimated S matrices.
    2. Calculate the mis-classification rate for each 30 points on the grid.  
    3. Among the grid points where the estimated L's rank is equal to 2, choose a pair of tuning parameters which gives the lowest mis-classification rate. 
```R
##### Model Selection through HNCV #####
gamma = seq(from=0.012,to=0.0128,by=0.0002);
delta = seq(from=0.04,to=0.05,by=0.002);

Kar_model <- Model_Sel(Karate_ad, gamma, delta);

count <- 1;
MisCl_rate1 <- matrix(rep(0,length(gamma)*length(delta)),nrow=length(gamma),ncol=length(delta));
for(i in 1:length(gamma)){
  for(j in 1:length(delta)){
    MisCl_rate1[i,j] <- CV(Karate_ad,gamma[i],delta[j],10)
    print(count)
    count <- count + 1;
  }
}

CV_ind <- which(MisCl_rate1==min(MisCl_rate1[Kar_model[[4]][1:5,]==2]),arr.ind=TRUE)
```
- Following lines of code evaluate the selected model through Network cross-validation.
Since the ground truth labels of nodes in the network are known, we can know true ad-hoc edges of the network.
  1. It is found that the selected model renders 11 pairs of ad-hoc edges. When compared with true ad-hoc edges, the model gives all 10 true ad-hoc edges with 1 extra edge (25-28).
  2. Six nodes are mis-classified as the Second community. (5th, 6h, 7th, 11th, 12th and 17th) 

```R
#### Selected Model
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

##### 11 pairs of ad-hoc edges are discovered #####
##### (1,9), (1,32), (2,31), (3,9), (3,10), (3,28), (3,29), (3,33)
##### (14,34), (20,34), (25,28)

##### Cluster nodes through K-means clustering #####
K=qr(Kara[[3]])$rank
vc <- eigen(Kara[[3]])$vectors[,1:2]
plot(vc[,1],vc[,2])

a<-kmeans(vc[,1:2], 2, iter.max = 1000, nstart = 100, 
       algorithm = "Hartigan-Wong")$cluster

C1 <- which(a==1,arr.ind=TRUE) 
# Nodes selected for C1 : 1  2  3  4  8 13 14 18 20 22
C2 <- which(a==2,arr.ind=TRUE) 
# Nodes selected for C2 : 5  6  7  9 10 11 12 15 16 17 19 21 23 24 25 26 27 28 29 30 31 32 33 34
```
- Here we draw the plot of network. Nodes classified as first community are colored in blue, nodes classified as second community are colored in orange. Discovered ad-hoc edges are colored in green. 
![](Codes/Karate.png)
