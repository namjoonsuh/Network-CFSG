setwd('C:/Users/nsuh3/Desktop/Citation Network/FLag/Codes & Data')
rm(list=ls())
require(igraph) # version 0.7.1. Note there are many changes in version 1.0.0
source('functions.R') #set the working directory

# read the data
authorPaperBiadj = as.matrix(read.table(file="authorPaperBiadj.txt",sep="\t", header=F))
authorList = as.matrix(read.table(file="authorList.txt",sep="\t", header=F, stringsAsFactors=F))
paperCitAdj = as.matrix(read.table("paperCitAdj.txt", header=F))
paperList = read.table("paperList.txt", sep=",", stringsAsFactors=F, header=T)

# compute adjacency matrix of the coauthorship network
coauthorAdjWeighted = authorPaperBiadj %*% t(authorPaperBiadj) 
coauthorAdj = (coauthorAdjWeighted >= 1) -0 
diag(coauthorAdj) = 0  

# construct the network (A)
thresh = 2  # the minimum number of coauthored papers for two authors to have an edge
# compute the adjacency matrix of the network (A)
coauthorAdjThresh = (coauthorAdjWeighted >= thresh) 
diag(coauthorAdjThresh) = 0
# find and label all components
ix = clusters(graph.adjacency(coauthorAdjThresh))
componentLabel = ix$membership
# find the labels for the largest 8 components
order(ix$csize, decreasing=T)[1:12] #output:  [1]   38  113  309  278  382  223 1036   42  175  725  996  117
# find their sizes
sort(ix$csize, decreasing=T)[1:12] #output:  [1] 236  18  14  13  10   9   9   8   8   8   8   7
#write.table(componentLabel, "data/coauthorThresh2ComponentLabel.txt", col.names=F, row.names=F)

# analyze the giant component (labeled as 38)
adj = coauthorAdjThresh[componentLabel == 38, componentLabel == 38]
author = authorList[componentLabel == 38]
# list 15 high-degree nodes
cat(author[order(rowSums(adj), decreasing=T)[1:15]], sep=", ")
# save the adjacency matrix of the giant component 
#write.table(adj, file="coauthorThresh2GiantAdj.txt", col.names=F,row.names=F)
#write.table(author, file="coauthorThresh2GiantAuthorList.txt", col.names=F,row.names=F)

# Fig 5 left: scree plot of the giant component 
plotScree(adj)

# plot the giant component
N = ncol(adj)
# pdf("coauthor-graph-thresh-hi-dim2.pdf", width=7, height=3.5)
par(mar=c(0,0,0,0))
set.seed(0)
vertexLabel = author
vertexLabel[apply(adj, 1, "sum")<7] =""
plot(graph.adjacency(adj, mode="undirected"), 
     vertex.color="white",
     vertex.label= vertexLabel,   vertex.label.cex=1.2, 
     vertex.size=3, vertex.label.color='red', asp=.5,
     edge.color='green', asp=.5, margin=c(0,0,0,0))
#dev.off()

############################# Model selection ##########################
############## Function of ADMM algorithm for estimation ###############
setwd('C:/Users/nsuh3/Desktop/Citation Network/FLag')
source('ADMM_Optim.R')
N <- nrow(adj)
gamma = seq(from=0.002124,to=0.002125,by=0.0000001);
delta = seq(from=0.009,to=0.01,by=0.0001);
lambda = 1; Count = 1;
AIC <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
BIC <- matrix( 0, nrow=length(gamma),ncol=length(delta) )

non_zeroS <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
L_Rank <- matrix( 0, nrow=length(gamma),ncol=length(delta) )
Like <- matrix( 0, nrow=length(gamma),ncol=length(delta) )

for(g in 1:length(gamma)){
  for(d in 1:length(delta)) {
    ### Use the ADMM method to estimate the parameters ###
    result <- ADMM(adj, gamma[g], delta[d])
    
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
    log_max <- a*sum(adj[upper.tri(adj,diag=FALSE)]) +
      (1/2)*sum(adj*(L+S)) - sum(log_sum[upper.tri(log_sum,diag=FALSE)])
    
    AIC[g,d] <- -2*log_max + M_absolute*2
    BIC[g,d] <- -2*log_max + M_absolute*(log((N*(N-1))/2))
    Like[g,d] <- log_max
    non_zeroS[g,d] <- non_zero_S
    L_Rank[g,d] <- qr(L)$rank
    
    print(Count)
    Count <- Count + 1
  }
}
###########################################################
result <- ADMM(adj, gamma[6], delta[11])
a1 <- result[[1]]
M1 <- result[[2]]
L1 <- result[[3]]
S1 <- result[[4]]

K <- qr(L)$rank

non_zero_S <- 0
for(i in 1:N){
  for(j in 1:N){
    if(i<j){
      if(S1[i,j] != 0)
        non_zero_S = non_zero_S + 1
    }
  }
}

############################## Community detection ############################
################ Differentiate communities by eyes ############################
DF = data.frame(X1=eigen(L1)$vectors[,1], X2=eigen(L1)$vectors[,2], Community = rep(0,nrow(adj)))
Comm_1 = intersect(which(DF[,1] <= -0.01, arr.ind = TRUE),which(DF[,2] <= -0.01, arr.ind = TRUE))
Comm_2 = intersect(which(DF[,1] >= 0.01, arr.ind = TRUE),which(DF[,2] <= -0.01, arr.ind = TRUE))

North_Carolina = author[Comm_1] # 25 researchers
Carrroll_Hall = author[Comm_2] # 37 researchers
Overlap_Commu = author[-c(Comm_1,Comm_2)] # 174 researchers

DF[Comm_1,3]=1; DF[Comm_2,3]=2; DF[-c(Comm_1,Comm_2),3]=3;
KMeans = kmeans(eigen(L1)$vectors[,1:K], K, iter.max = 1000, nstart = 100, algorithm = "Hartigan-Wong")
plot(DF[,1:2],col=DF$Community,pch=2,xlab="First Eigen",ylab="Second Eigen")
legend(-0.1, 0.1, legend=c("Carroll-Hall", "Mixed Community","North Carolina"),
       col = c("red","green","black"), pch=2, cex=0.8)

author = cbind(author,rep(0,length(author)))
author[Comm_1,2]="North Carolina"
author[Comm_2,2]="Carroll-Hall"
author[-c(Comm_1,Comm_2),2]="Mixed Community"

################ Differentiate communities by K-means ##############
K = 2
KMeans = kmeans(eigen(L1)$vectors[,1:K], K, iter.max = 1000, nstart = 100, algorithm = "Hartigan-Wong")
C1 = which(KMeans$cluster == 1, arr.ind = TRUE)
C2 = which(KMeans$cluster == 2, arr.ind = TRUE)

############################## Ad-hoc dependency ############################
Edge = which(S1>0, arr.ind=TRUE)
g = graph_from_edgelist(Edge)
S_g = as.undirected(g)
Edge_list = get.edgelist(S_g)

Ad_hoc_link <- data.frame(Pair1 = Edge_list[,1], Pair2 = Edge_list[,2],  S_val = S1[Edge_list])
Sorted_Ad_hoc_link <- Ad_hoc_link[order(-Ad_hoc_link$S_val),] 
Edge_list <- as.matrix(Sorted_Ad_hoc_link[,1:2])

row_name = row.names(adj)
col1 = Edge_list[,1]
col2 = Edge_list[,2]
C1 = author[col1,]
C2 = author[col2,]
cbind(C1,C2)



















set.seed(1233)
K = 2
KMeans = kmeans(eigen(L1)$vectors[,1:K], K, iter.max = 1000, nstart = 100, algorithm = "Hartigan-Wong")
Topic1 = which(KMeans$cluster == 1, arr.ind = TRUE)
Topic2 = which(KMeans$cluster == 2, arr.ind = TRUE)
Topic3 = which(KMeans$cluster == 3, arr.ind = TRUE)

paperList = read.table("C:/Users/nsuh3/Desktop/Citation Network/paperList.txt", sep=",", stringsAsFactors=F, header=T)
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

K=3
DF <- data.frame(FE = eigen(L1)$vector[,1], SE = eigen(L1)$vector[,2], Cluster = rep(0,N))
DF[Topic1,3]=1; DF[Topic2,3]=2; DF[Topic3,3]=3; 
KMeans = kmeans(eigen(L1)$vectors[,1:K], K, iter.max = 1000, nstart = 100, algorithm = "Hartigan-Wong")
plot(DF[,1:2], col=KMeans$cluster,pch=2,xlab="First Eigen",ylab="Second Eigen")
legend(-0.3, 0.33, legend=c("Multiple testing", "Variable Selection", "Mixed Topic"),
       col = c("red","green","black"), pch=2, cex=0.8)


Edge = which(S>0, arr.ind=TRUE)
g = graph_from_edgelist(Edge)
S_g = as.undirected(g)

Ad_hoc_S <- graph_from_edgelist(el, directed = TRUE)
plot(Ad_hoc_S)

A=get.edgelist(adj)
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

E(X_draw)$color <- "black"
E(X_draw)[list]$color <- "blue"
plot(X_draw)
