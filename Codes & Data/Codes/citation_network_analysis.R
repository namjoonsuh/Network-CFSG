####################  preparation  ####################
rm(list=ls())
setwd('C:/Users/namjo/Documents/GitHub/Citation-Network/Codes & Data/Codes')
library('Matrix')      ## For matrix computation on BIC
library('igraph')      ## For drawing graph
source('SynData.R')    ## Function for creating binary network
source('ADMM_Optim.R') ## Function of ADMM algorithm for estimation
source('functions.R')
source('GD.R')
################### Load the Citation Network Data ###################
setwd('C:/Users/namjo/Documents/GitHub/Citation-Network/Codes & Data/Statistician network data')
paperCitAdj = as.matrix(read.table("paperCitAdj.txt", header=FALSE))
G <- graph.adjacency(paperCitAdj, mode="undirected", weighted=NULL) # For directed networks
X <- as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE,attr=NULL)
Node_deg <- rep(0,nrow(X))
for (i in 1:sqrt(length(X))){
  Node_deg[i]<-sum(X[i,])
}
###### Get a network whose degree is greater than or equal to 7 ######
ind <- which(Node_deg>=7,arr.ind=TRUE)
X_new <- X[ind,ind]
N <- ncol(X_new)

X_draw <- graph_from_adjacency_matrix(X_new, mode = c("undirected"))
plot(X_draw)
plotScree(X_new)
X = X_new

################### Model Selection ###################
gamma = seq(from=0.01,to=0.01,by=0.00000);
delta = seq(from=0.01,to=0.01,by=0.0001);
lambda = 1; Count = 1;
AIC <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
BIC <- matrix( 0, nrow=length(gamma),ncol=length(delta) )

non_zeroS <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
L_Rank <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
Like <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
N = ncol(X_new)

source('ADMM_Optim.R') ## Function of ADMM algorithm for estimation
for(g in 1:length(gamma)){
  for(d in 1:length(delta)) {
    ### Use the ADMM method to estimate the parameters ###
    result <- ADMM(X, 1, 0.01)
    
    a<-result1[[1]]
    M<-result1[[2]]
    L<-result1[[3]]
    S<-result1[[4]]
    
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
    log_max <- a*sum(X[upper.tri(X,diag=FALSE)]) +
      (1/2)*sum(X*(L+S)) - sum(log_sum[upper.tri(log_sum,diag=FALSE)])
    
    AIC[g,d] <- -2*log_max + M_absolute*2
    BIC[g,d] <- -2*log_max + M_absolute*(log((N*(N-1))/2))
    Like[g,d] <- log_max
    non_zeroS[g,d] <- non_zero_S
    L_Rank[g,d] <- qr(L)$rank
    
    print(Count)
    Count <- Count + 1
  }
}

#################### Choose the proper tuning parameters #######################
result1 <- ADMM(X, 0.000912, 0.0097)  ## Gamma : 0.000912, Delta : 0.0097
a1<-result1[[1]]
M1<-result1[[2]]
L1<-result1[[3]]
S1<-result1[[4]]

######### Calculate the rank of L and count the number of non-zero elements in S #########
K = qr(L1)$rank
N = ncol(X)

non_zero_S <- 0
for(i in 1:N){
  for(j in 1:N){
    if(i<j){
      if(S1[i,j] != 0)
        non_zero_S = non_zero_S + 1
    }
  }
}

####################### Cluster the nodes by K-means #######################
set.seed(1233)
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

####################### Further Analysis on T3 #######################
############# Analysis on another big chunk of data ##################
source('functions.R')
par(mfrow=c(1,2),mar = c(2, 2, 2, 1))
plotScree(X)
X_chunk = X[Topic3,Topic3]
plotScree(X_chunk)
# I guess.. there are 3 topics and also another chunk of papers with mixed topics
############# Model Selection #############
gamma = seq(from=0.0015,to=0.002,by=0.0001);
delta = seq(from=0.0062,to=0.008,by=0.0002);
lambda = 1; Count = 1;
AIC <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
BIC <- matrix( 0, nrow=length(gamma),ncol=length(delta) )

non_zeroS <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
L_Rank <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
Like <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
N = ncol(X_chunk)

source('ADMM_Optim.R') ## Function of ADMM algorithm for estimation
for(g in 1:length(gamma)){
  for(d in 1:length(delta)) {
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
    M_absolute <- (N*K-((K-1)*K/2)) + non_zero_S + 1
    
    log_sum<-log(1+exp(a + L + S))
    
    log_max <- 0;
    log_max <- a*sum(X_chunk[upper.tri(X_chunk,diag=FALSE)]) +
      (1/2)*sum(X_chunk*(L+S)) - sum(log_sum[upper.tri(log_sum,diag=FALSE)])
    
    AIC[g,d] <- -2*log_max + M_absolute*2
    BIC[g,d] <- -2*log_max + M_absolute*(log((N*(N-1))/2))
    Like[g,d] <- log_max
    non_zeroS[g,d] <- non_zero_S
    L_Rank[g,d] <- qr(L)$rank
    
    print(Count)
    Count <- Count + 1
  }
}
############################# Choose the proper tuning parameters #############################
result2 <- ADMM(X_chunk, gamma[3], delta[8]) 
N = ncol(X_chunk)
a2<-result2[[1]]
M2<-result2[[2]]
L2<-result2[[3]]
S2<-result2[[4]]
############# Clusters by Kmeans #############
K=4
set.seed(1232)
KMeans = kmeans(eigen(L2)$vectors[,1:K], K, iter.max = 1000, nstart = 100, algorithm = "Hartigan-Wong")
Topic3 = which(KMeans$cluster == 1,arr.ind = TRUE)
Topic4 = which(KMeans$cluster == 2,arr.ind = TRUE)
Topic5 = which(KMeans$cluster == 3,arr.ind = TRUE)
Topic6 = which(KMeans$cluster == 4,arr.ind = TRUE)

paperList = read.table("C:/Users/nsuh3/Desktop/Citation Network/paperList.txt", sep=",", stringsAsFactors=F, header=T)
community = rep(0,nrow(paperList))
paperList = cbind(paperList,community)

row_name = row.names(X_chunk)
Topic_3 = as.numeric(sub('.','',row_name[Topic3]))
Topic_4 = as.numeric(sub('.','',row_name[Topic4]))
Topic_5 = as.numeric(sub('.','',row_name[Topic5]))
Topic_6 = as.numeric(sub('.','',row_name[Topic6]))

pap_list = paperList[,3]
T3 = pap_list[Topic_3]
T4 = pap_list[Topic_4]
T5 = pap_list[Topic_5]
T6 = pap_list[Topic_6]

N = ncol(X_chunk)
DF <- data.frame(FE = eigen(L2)$vector[,1], SE = eigen(L2)$vector[,2], Cluster = rep(0,N))
DF[Topic3,3]=1; DF[Topic4,3]=2; DF[Topic5,3]=3; DF[Topic6,4]=4; 
plot(DF[,1:2], col=KMeans$cluster,pch=2,xlab="First Eigen",ylab="Second Eigen")
legend(-0.05, 0.28, legend=c("Functional Analysis", "Dimension Reduction", "Bayesian Statistics", "Mixed Topic"),
       col = c("black","red","green","blue"), pch=2, cex=0.8)

paperList[Topic_3,6] = "Bayesian Statistics"
paperList[Topic_4,6] = "Functional Analysis"
paperList[Topic_5,6] = "Mixed Topics"
paperList[Topic_6,6] = "Dimension Reduction"

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
write.csv(Ad_hoc_list, "C:/Users/nsuh3/Desktop/Citation Network/FLag/ad_hoc_192.csv") 

################################################################################
############################ Get the list of Topics ############################
X_Topic1 = X[Topic1,Topic1]
X_draw_Topic1 <- graph_from_adjacency_matrix(X_Topic1, mode = c("undirected"))
plot(X_draw_Topic1, main = "Papers on Multiple Hypothesis Testing")

X_Topic2 = X[Topic3,Topic3]
X_draw_Topic2 <- graph_from_adjacency_matrix(X_Topic2, mode = c("undirected"))
plot(X_draw_Topic2, main = "Papers on variable selection")

######## Hui Zou's "On the degrees of freedom of LASSO" ##########
H_Zou = which(row.names(X_new)=="V507",arr.ind=TRUE)
H_Zou_Cit_pap = which(X_new[H_Zou,]==1, arr.ind=TRUE)
H_Zou_Cit_pap_ind = as.numeric(sub('.','',row.names(X_new)[H_Zou_Cit_pap]))
H_Zou_paper = paperList[H_Zou_Cit_pap_ind,c(3,6)]
length(which(H_Zou_paper[,2]=="Mixed Topics",arr.ind=TRUE))
##################################################################
