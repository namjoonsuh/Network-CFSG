####################  preparation  ####################
rm(list=ls())
library('Matrix')      ## For matrix computation on BIC
library('igraph')      ## For drawing graph

################### Load the Citation Network Data ###################
setwd('C:/Users/Namjoon Suh/Desktop/Stat+PDE/Citation-Network/Codes & Data/Statistician network data')
paperCitAdj = as.matrix(read.table("paperCitAdj.txt", header=FALSE))
G <- graph.adjacency(paperCitAdj, mode="undirected", weighted=NULL) # For directed networks
X <- as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE,attr=NULL)
Node_deg <- rep(0,nrow(X))
for (i in 1:sqrt(length(X))){
  Node_deg[i]<-sum(X[i,])
}

################### Load Functions for analysis ###################
setwd('C:/Users/Namjoon Suh/Desktop/Stat+PDE/Citation-Network/Codes & Data/Codes')
source('SynData.R')    ## Function for creating binary network
source('ADMM_Optim.R') ## Function of ADMM algorithm for estimation
source('functions.R')
source('GD.R')

###### Get a network whose degree is greater than or equal to 10 ######
ind <- which(Node_deg>=10,arr.ind=TRUE)
X_new <- X[ind,ind]
N <- ncol(X_new)
A = rowSums(X_new)
ind = which(A==0,arr.ind=TRUE)
X = X_new[-ind,-ind]

par(mfrow=c(1,2),mar = c(2, 2, 1.5, 1.5)) # ??? ??? ??? ??? # 
plotScree(X)
title(TeX("$\\X^{orig}$"),cex.main=0.85)

plotScree(X_chunk)
title(TeX("$\\X^{sub}$"),cex.main=0.85)

################### Model Selection ###################
gamma = seq(from=0.0021091,to=0.0021095,by=0.00000006);
delta = seq(from=0.01913,to=0.0194,by=0.00003);
lambda = 1; Count = 1;
AIC <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
BIC <- matrix( 0, nrow=length(gamma),ncol=length(delta) )

non_zeroS <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
L_Rank <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
Like <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
N = ncol(X)

setwd('C:/Users/namjo/Documents/GitHub/Citation-Network/Codes & Data/Codes') ## Function of ADMM algorithm for estimation
source('ADMM_Optim.R') ## Function of ADMM algorithm for estimation
source('functions.R')
source('GD.R')

for(g in 1:length(gamma)){
  for(d in 1:length(delta)) {
    ### Use the ADMM method to estimate the parameters ###
    result <- ADMM(X, gamma[g], delta[d])

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
    
    non_zeroS[g,d] <- non_zero_S
    L_Rank[g,d] <- qr(L)$rank
    
    print(Count)
    Count <- Count + 1
  }
}
#################### Choose the proper tuning parameters #######################
N = ncol(X)
result1 <- ADMM(X, 0.002109, 0.0194)  ## Gamma : 0.002109, Delta : 0.0194
a1 <- result1[[1]]
M1 <- result1[[2]]
L1 <- result1[[3]]
S1 <- result1[[4]]

K<-qr(L1)$rank; K;
non_zero_S <- 0
for(i in 1:N){
  for(j in 1:N){
    if(i<j){
      if(S1[i,j] != 0)
        non_zero_S = non_zero_S + 1
    }
  }
}
non_zero_S;

####################### Cluster the nodes by K-means #######################
set.seed(1233)
K = qr(L1)$rank;
setwd('C:/Users/namjo/Documents/GitHub/Citation-Network/Codes & Data/Statistician network data')
KMeans = kmeans(eigen(L1)$vectors[,1:K], K, iter.max = 1000, nstart = 100, algorithm = "Hartigan-Wong")
Topic1 = which(KMeans$cluster == 1,arr.ind = TRUE)
Topic2 = which(KMeans$cluster == 2,arr.ind = TRUE)
Topic3 = which(KMeans$cluster == 3,arr.ind = TRUE)

paperList = read.table("paperList.txt", sep=",", stringsAsFactors=F, header=T)
community = rep(0,nrow(paperList))
paperList = cbind(paperList,community)

row_name = row.names(X)
Topic_1 = as.numeric(sub('.','',row_name[Topic1]))
Topic_2 = as.numeric(sub('.','',row_name[Topic2]))
Topic_3 = as.numeric(sub('.','',row_name[Topic3]))

pap_list = paperList[,3]
T1 = pap_list[Topic_1]
T2 = pap_list[Topic_2]
T3 = pap_list[Topic_3]

paperList[Topic_1,6] = "Variable Selection"
paperList[Topic_2,6] = "Multiple Testing"

################################################################################
############################ Get the list of Topics ############################
################################################################################
X_Topic1 = X[Topic1,Topic1]
X_draw_Topic1 <- graph_from_adjacency_matrix(X_Topic1, mode = c("undirected"))
plot(X_draw_Topic1, main = "Papers on Multiple Hypothesis Testing")

X_Topic2 = X[Topic2,Topic2]
X_draw_Topic2 <- graph_from_adjacency_matrix(X_Topic2, mode = c("undirected"))
plot(X_draw_Topic2, main = "Papers on variable selection")

X_Topic3 = X[Topic3,Topic3]
X_draw_Topic3 <- graph_from_adjacency_matrix(X_Topic3, mode = c("undirected"))
plot(X_draw_Topic3, main = "Papers on variable selection")

sum(X)/2
999-(sum(X_Topic1)/2 + sum(X_Topic2)/2 + sum(X_Topic3)/2)
non

X_Topic4 = X[Topic4,Topic4]
X_draw_Topic4 <- graph_from_adjacency_matrix(X_Topic4, mode = c("undirected"))
plot(X_draw_Topic4, main = "Papers on variable selection")
par(mfrow=c(1,2),mar = c(4, 4, 2, 1))

DF <- data.frame(FE = eigen(L1)$vector[,1], SE = eigen(L1)$vector[,2], Cluster = rep(0,N))
DF[Topic1,3]=1; DF[Topic2,3]=2; DF[Topic3,3]=3;
plot(DF[,1:2], col=KMeans$cluster,pch=2,xlab="First Eigen",ylab="Second Eigen")
legend(-0.3, 0.33, legend=c("Variable Selection", "Mixed Clusters", "Multiple Testing"),
       col = c("red","green","black"), pch=2, cex=0.8)

####################### See the lists of papers for each topics #######################
T1 ### Variable Selection
T2 ### Multiple Testing (FDR control)
T3 ### Mixed Topic
T4 ### 

########################################################################################
############################ Ad-hoc dependent structure : S ############################
################## Getting the paper list & Analysis on ad-hoc link structures ###############
Edge = which(S1>0, arr.ind=TRUE)
g = graph_from_edgelist(Edge)
S_g = as.undirected(g)
Edge_list = get.edgelist(S_g)

Ad_hoc_link <- data.frame(Pair1 = Edge_list[,1], Pair2 = Edge_list[,2],  S_val = S1[Edge_list])
Sorted_Ad_hoc_link <- Ad_hoc_link[order(-Ad_hoc_link$S_val),] 
Edge_list <- as.matrix(Sorted_Ad_hoc_link[,1:2])

row_name = row.names(X)
col1 = as.numeric(sub('.','',row_name[Edge_list[,1]]))
col2 = as.numeric(sub('.','',row_name[Edge_list[,2]]))
Paper_Table = matrix(as.numeric(cbind(col1,col2)),length(Edge_list)/2,2)
pap_list = paperList[,3]
C1 = pap_list[col1]
C2 = pap_list[col2]

Ad_hoc_list = list()
P1 = paperList[col1,c(1,6,3)]
P2 = paperList[col2,c(1,6,3)]
for(i in 1:length(col1)){
  Ad_hoc_list = rbind(Ad_hoc_list,P1[i,],P2[i,])
}
Ad_hoc_list = cbind(10^4*Sorted_Ad_hoc_link[,3],Ad_hoc_list)
Ad_hoc_list
################# View the papers of Ad-hoc pair #################
which(Ad_hoc_list[,3]=="Multiple Testing")
which(Ad_hoc_list[,3]=="Variable Selection")
write.csv(Ad_hoc_list, "C:/Users/namjo/Documents/GitHub/Citation-Network/Codes & Data/Codes/ad_hoc_192.csv") 

X_draw <- graph_from_adjacency_matrix(X, mode = c("undirected"))
V(X_draw)[Topic1]$color <- "orange"
V(X_draw)[Topic2]$color <- "black"
V(X_draw)[Topic3]$color <- "grey"
V(X_draw)[Topic4]$color <- "blue"
V(X_draw)[Topic5]$color <- "purple"
V(X_draw)[Topic6]$color <- "red"
V(X_draw)[Topic7]$color <- "green"
V(X_draw)[Topic8]$color <- "yellow"

E(X_draw)$color <- "red"
E(X_draw)[Edge_list]$color <- "green"
plot(X_draw, vertex.size=10, vertex.label=NA,layout=layout_with_fr)

par(mfrow=c(1,1))

####################### Further Analysis on T3 #######################
############# Analysis on another big chunk of data ##################
source('functions.R')
par(mfrow=c(1,2),mar = c(2, 2, 2, 1))
plotScree(X)
X_chunk = X[Topic3,Topic3]
plotScree(X_chunk)
title(TeX("$\\X^{sub}$"),cex.main=1)
# I guess.. there are 3 topics and also another chunk of papers with mixed topics
############# Model Selection #############
gamma = seq(from=0.0031,to=0.00315,by=0.000005);
delta = seq(from=0.014,to=0.015,by=0.0001);
lambda = 1; Count = 1;
AIC <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
BIC <- matrix( 0, nrow=length(gamma),ncol=length(delta) )

non_zeroS <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
L_Rank <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
Like <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
N = ncol(X_chunk)

setwd('C:/Users/namjo/Documents/GitHub/Citation-Network/Codes & Data/Codes') ## Function of ADMM algorithm for estimation
source('ADMM_Optim.R') ## Function of ADMM algorithm for estimation
source('functions.R')
source('GD.R')
for(g in 1 : length(gamma)){
  for(d in 1 : length(delta)) {
    ### Use the ADMM method to estimate the parameters ###
    result <- ADMM(X_chunk, gamma[g], delta[d])
    
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
    
    non_zeroS[g,d] <- non_zero_S
    L_Rank[g,d] <- qr(L)$rank
    
    print(Count)
    Count <- Count + 1
  }
}

############################# Choose the proper tuning parameters #############################
result2 <- ADMM(X_chunk, 0.00312, 0.0146) # Gamma : 0.00312, delta : 0.0146
N = ncol(X_chunk)
a2<-result2[[1]]
M2<-result2[[2]]
L2<-result2[[3]]
S2<-result2[[4]]

############# Clusters by Kmeans #############
K=qr(L2)$rank;
KMeans = kmeans(eigen(L2)$vectors[,1:K], K, iter.max = 1000, nstart = 100, algorithm = "Hartigan-Wong")
Topic3 = which(KMeans$cluster == 1,arr.ind = TRUE)
Topic4 = which(KMeans$cluster == 2,arr.ind = TRUE)
Topic5 = which(KMeans$cluster == 3,arr.ind = TRUE)
Topic6 = which(KMeans$cluster == 4,arr.ind = TRUE)
Topic7 = which(KMeans$cluster == 5,arr.ind = TRUE)

paperList = read.table("C:/Users/nsuh3/Desktop/Citation Network/paperList.txt", sep=",", stringsAsFactors=F, header=T)
community = rep(0,nrow(paperList))
paperList = cbind(paperList,community)

row_name = row.names(X_chunk)
Topic_3 = as.numeric(sub('.','',row_name[Topic3]))
Topic_4 = as.numeric(sub('.','',row_name[Topic4]))
Topic_5 = as.numeric(sub('.','',row_name[Topic5]))
Topic_6 = as.numeric(sub('.','',row_name[Topic6]))
Topic_7 = as.numeric(sub('.','',row_name[Topic7]))

pap_list = paperList[,3]
T3 = pap_list[Topic_3]
T4 = pap_list[Topic_4]
T5 = pap_list[Topic_5]
T6 = pap_list[Topic_6]
T7 = pap_list[Topic_7]

N = ncol(X_chunk)
DF <- data.frame(FE = eigen(L2)$vector[,1], SE = eigen(L2)$vector[,2], Cluster = rep(0,N))
DF[Topic3,3]=1; DF[Topic4,3]=2; DF[Topic5,3]=3; DF[Topic6,4]=4; 
plot(DF[,1:2], col=KMeans$cluster,pch=2,xlab="First Eigen",ylab="Second Eigen")
legend(-0.05, 0.28, legend=c("Functional Analysis", "Dimension Reduction", "Bayesian Statistics", "Mixed Topic"),
       col = c("black","red","green","blue"), pch=2, cex=0.8)

paperList[Topic_3,6] = "Mixed Topics"
paperList[Topic_4,6] = "Dimension Reduction"
paperList[Topic_5,6] = "Bayesian Nonparametric Stat"
paperList[Topic_6,6] = "Covariance Estimation"
paperList[Topic_7,6] = "Functional/Longitudinal Analysis"

########################################################################################
############################ Ad-hoc dependent structure : S ############################
############# Getting the paper list & Analysis on ad-hoc link structures ##############
Edge = which(S1>0, arr.ind=TRUE)
g = graph_from_edgelist(Edge)
S_g = as.undirected(g)
Edge_list = get.edgelist(S_g)

Ad_hoc_link <- data.frame(Pair1 = Edge_list[,1], Pair2 = Edge_list[,2],  S_val = S1[Edge_list])
Sorted_Ad_hoc_link <- Ad_hoc_link[order(-Ad_hoc_link$S_val),] 
Edge_list <- as.matrix(Sorted_Ad_hoc_link[,1:2])

row_name = row.names(X)
col1 = as.numeric(sub('.','',row_name[Edge_list[,1]]))
col2 = as.numeric(sub('.','',row_name[Edge_list[,2]]))
Paper_Table = matrix(as.numeric(cbind(col1,col2)),length(Edge_list)/2,2)
pap_list = paperList[,3]
C1 = pap_list[col1]
C2 = pap_list[col2]

Ad_hoc_list = list()
P1 = paperList[col1,c(1,6,3)]
P2 = paperList[col2,c(1,6,3)]
for(i in 1:length(col1)){
  Ad_hoc_list = rbind(Ad_hoc_list,P1[i,],P2[i,])
}
Ad_hoc_list = cbind(10^4*Sorted_Ad_hoc_link[,3],Ad_hoc_list)
Ad_hoc_list
################# View the papers of Ad-hoc pair #################
which(Ad_hoc_list[,3]=="Multiple Testing")
which(Ad_hoc_list[,3]=="Variable Selection")
write.csv(Ad_hoc_list, "C:/Users/namjo/Documents/GitHub/Citation-Network/Codes & Data/Codes/ad_hoc_51.csv") 

################################################################################
############################ Get the list of Topics ############################
################################################################################
X_Topic1 = X[Topic1,Topic1]
X_draw_Topic1 <- graph_from_adjacency_matrix(X_Topic1, mode = c("undirected"))
plot(X_draw_Topic1, main = "Papers on Multiple Hypothesis Testing")

X_Topic2 = X[Topic2,Topic2]
X_draw_Topic2 <- graph_from_adjacency_matrix(X_Topic2, mode = c("undirected"))
plot(X_draw_Topic2, main = "Papers on variable selection")

X_Topic3 = X[Topic3,Topic3]
X_draw_Topic3 <- graph_from_adjacency_matrix(X_Topic3, mode = c("undirected"))
plot(X_draw_Topic3, main = "Papers on variable selection")

X_Topic4 = X[Topic4,Topic4]
X_draw_Topic4 <- graph_from_adjacency_matrix(X_Topic4, mode = c("undirected"))
plot(X_draw_Topic4, main = "Papers on variable selection")

######## Hui Zou's "On the degrees of freedom of LASSO" ##########
H_Zou = which(row.names(X)=="V507",arr.ind=TRUE)
H_Zou_Cit_pap = which(X[H_Zou,]==1, arr.ind=TRUE)
H_Zou_Cit_pap_ind = as.numeric(sub('.','',row.names(X)[H_Zou_Cit_pap]))
H_Zou_paper = paperList[H_Zou_Cit_pap_ind,c(3,6)]
length(which(H_Zou_paper[,2]=="Mixed Topics",arr.ind=TRUE))
##################################################################
