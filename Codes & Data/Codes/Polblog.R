rm(list=ls())
setwd('C:/Users/namjo/Desktop/Research Files/Network Paper/Codes & Data/Codes')
source('ADMM_Optim.R') # We set r = 4 in ADMM function, delta = 10^-7
source('Synthetic Networks.R')
library("igraph"); library("CINNA");

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

## Algorithm of Ma et al paper, 2020.
Init<-function(A,X,eta,K,gamma,Iter){
  n <- nrow(A); ones <- t(t(rep(1,n))); 
  G <- matrix(0,n,n); al <- t(t(rep(0,n)));
  lambda <- 2*sqrt(sum(A)/n)
  
  for(t in 1:Iter){
    Theta_t <- al%*%t(ones) + ones%*%t(al) + G;
    sig_T = 1/(1+exp(-Theta_t))
    G = G + eta * (A - sig_T - lambda*diag(n) - gamma*G);
    al = al + eta * ((A - sig_T)%*%ones/(2*n) - gamma*al);
    #b = b + eta * ((A - sig_T)/sum(X^2)  - gamma*b) 
  }
  
  E_Val <- eigen(G)$values[1:K]; E_Vec <- eigen(G)$vectors[,1:K];
  Z = E_Vec%*%diag(sqrt(E_Val))
  res <- list(Z,al)
  return(res)
}
Algo<-function(A,X,K,res,eta, Iter){
  n <- nrow(A); ones <- t(t(rep(1,n))); 
  Z = res[[1]]; al = res[[2]];
  eta_Z = eta/norm(Z,"O")^2; eta_al = eta/(2*n);
  
  for(t in 1:Iter){
    Theta_t <- al%*%t(ones) + ones%*%t(al) + Z%*%t(Z);
    sig_T = 1/(1+exp(-Theta_t));
    Z = Z + eta_Z*(A-sig_T)%*%Z;
    al = al + eta_al*(A-sig_T)%*%ones;
  }
  Res <- list(al,Z)
  
  return(Res)
}

## Load the Political Blog Data set ##
pol = as.matrix(read.table("polblogsgiant.txt"));
y = read.table("polblogsgianty.txt")
gamma = seq(from=0.00061,to=0.0007,by=0.00001);
delta = seq(from=0.00651,to=0.0075,by=0.0001);
Res_pol <- Model_Sel(pol,gamma,delta);
save(Res_pol, file = "PolBlog_result.RData")
load("PolBlog_result.RData")

### Run CV
count <- 1;
MisCl_rate_pol <- matrix(rep(0,length(gamma)^2),nrow=length(gamma),ncol=length(delta));
for(i in 1:10){
  for(j in 1:10){
    MisCl_rate_pol[i,j] <- CV(pol,gamma[i],delta[j],10)    
    print(count)
    count <- count + 1;
  }
}
save(MisCl_rate_pol, file = "MisCl_rate_pol_result.RData")
load("MisCl_rate_pol_result.RData")

### Pairs of parameters chosen by Network Cross-validation : gamma : 0.00063, delta : 0.00661
Res <- ADMM(pol,0.00058,0.0061)

non_zero_S <- 0
for(i in 1:1222){
  for(j in 1:1222){
    if(i<j){
      if(Res[[4]][i,j] > 0){
        non_zero_S = non_zero_S + 1
      }
    }
  }
}
non_zero_S

################## Experiment on Goodness of fit on data ##################
### 1. Fit with combined latent-logistic regression model ###
B=1000; Bootstrap <- rep(0,B);
for(i in 1:B){
  X = CrtData(Res[[1]],Res[[3]],Res[[4]],1222)[[1]]
  log_sum <- log(1+exp(Res[[1]] + Res[[3]] + Res[[4]]))
  
  log_max <- 0;
  log_max <- Res[[1]]*sum(X[upper.tri(X,diag=FALSE)]) +
    (1/2)*sum(X*(Res[[3]]+Res[[4]])) - sum(log_sum[upper.tri(log_sum,diag=FALSE)])
  Bootstrap[i] = log_max;
  print(i)
}

# pol_likelihood : -80012.56
log_sum <- log(1+exp(Res[[1]] + Res[[3]] + Res[[4]]))
pol_likelihood <- Res[[1]]*sum(pol[upper.tri(pol,diag=FALSE)]) +
  (1/2)*sum(pol*(Res[[3]]+Res[[4]])) - sum(log_sum[upper.tri(log_sum,diag=FALSE)])

hist(Bootstrap,breaks=100)
save(Bootstrap, file = "Combined_Latent_Logistics_Bootstrap.RData")

### 2. Fit with latent variable model ###
Res_1 <- ADMM(pol,1,0.0061)
qr(Res_1[[3]])$rank
non_zero_S <- 0
for(i in 1:1222){
  for(j in 1:1222){
    if(i<j){
      if(Res_1[[4]][i,j] != 0)
        non_zero_S = non_zero_S + 1
    }
  }
}

B=1000; Bootstrap_1 <- rep(0,B);
for(i in 1:B){
  X = CrtData(Res_1[[1]],Res_1[[3]],Res_1[[4]],1222)[[1]]
  log_sum <- log(1+exp(Res_1[[1]] + Res_1[[3]] + Res_1[[4]]))
  
  log_max <- 0;
  log_max <- Res_1[[1]]*sum(X[upper.tri(X,diag=FALSE)]) +
    (1/2)*sum(X*(Res_1[[3]]+Res_1[[4]])) - sum(log_sum[upper.tri(log_sum,diag=FALSE)])
  Bootstrap_1[i] = log_max;
  print(i)
}
hist(Bootstrap_1,breaks=100)
save(Bootstrap_1, file = "Latent_var_Bootstrap.RData")

# pol_likelihood_1 : -80012.57
log_sum <- log(1+exp(Res_1[[1]] + Res_1[[3]] + Res_1[[4]]))
pol_likelihood_1 <- Res_1[[1]]*sum(pol[upper.tri(pol,diag=FALSE)]) +
  (1/2)*sum(pol*(Res_1[[3]]+Res_1[[4]])) - sum(log_sum[upper.tri(log_sum,diag=FALSE)])

### 3. Fit with latent variable model in Ma et al ###
CrtData_Ma <- function(P,N){
  # Create empty matrix for storing adjacency matrix X and Probability for 
  X <- matrix(0,N,N); 

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
res <- Init(pol,0,0.2,2,0.001,10)
Result_Ma <- Algo(pol,0,2,res,0.2,10)
al <- Result_Ma[[1]]; Z <- Result_Ma[[2]];  ones <- t(t(rep(1,1222)));
P = exp(al%*%t(ones) + ones%*%t(al)+Z%*%t(Z))/(1+exp(al%*%t(ones) + ones%*%t(al)+Z%*%t(Z)))
log_sum_Ma <- log(1+exp(al%*%t(ones) + ones%*%t(al) + Z%*%t(Z)))

B=1000; Bootstrap_Ma <- rep(0,B);
for(i in 428:B){
  X = CrtData_Ma(P,1222)
  pol_likelihood_Ma <- (1/2)*sum(X*(al%*%t(ones) + ones%*%t(al) + Z%*%t(Z)))-
    sum(log_sum_Ma[upper.tri(log_sum_Ma,diag=FALSE)])
  Bootstrap_Ma[i] = pol_likelihood_Ma;
  print(i)
}
hist(Bootstrap_Ma,breaks=100)
pol_likelihood_Ma <- (1/2)*sum(pol*(al%*%t(ones) + ones%*%t(al) + Z%*%t(Z)))-
  sum(log_sum_Ma[upper.tri(log_sum_Ma,diag=FALSE)])
save(Bootstrap_Ma, file = "Latent_var_Bootstrap_Ma.RData")

load("Latent_var_Bootstrap.RData")
load("Latent_var_Bootstrap_Ma.RData")

mean <- mean(Bootstrap_1); sd <- sd(Bootstrap_1);
df <- data.frame(likelihood = Bootstrap_1)
binwidth <- 100; n = 1000;

Z <- (-80012.57-mean)/sd
pnorm(Z, 0, 1)

p1<-ggplot(df, aes(x = Bootstrap_1, mean = mean, sd = sd, binwidth = binwidth),n = n) +
  theme_bw() +
  geom_histogram(binwidth = binwidth, 
                 colour = "white", fill = "cornflowerblue", size = 0.1) +
  geom_vline(xintercept = -80012.57, linetype="dashed", 
             color = "red", size=1.5) +
  xlab(TeX('likelihood of LFaS-LR')) 

  # stat_function(fun = function(x) dnorm(x, mean = mean, sd = sd) * n * binwidth,
  #             color = "darkred", size = 1)

p2<-ggplot(df, aes(x = Bootstrap_Ma, mean = mean, sd = sd, binwidth = binwidth),n = n) +
  theme_bw() +
  geom_histogram(binwidth = binwidth, 
                 colour = "white", fill = "cornflowerblue", size = 0.1) +
  xlab(TeX('likelihood of Inn-Pr')) 
  # stat_function(fun = function(x) dnorm(x, mean = mean, sd = sd) * n * binwidth,
  #             color = "darkred", size = 1)


library(cowplot)
plot_grid(p1,p2,labels="AUTO")


par(xpd=FALSE); library(extrafont); loadfonts(device = "win");
par(family = "LM Roman 10");
par(mfrow=c(1,2))
hist(Bootstrap_1,breaks=100,main="Bootstrap Distribution")
abline(v = pol_likelihood, col="red", lwd=3, lty=2)

hist(Bootstrap_Ma,breaks=100,main="Bootstrap Distribution")
abline(v = pol_likelihood_Ma, col="red", lwd=3, lty=2)

####### True Positive #######
source('ADMM_Optim.R') # We set r = 4 in ADMM function.  
Res <- ADMM(pol,0.00063,0.0061)

###### |S*| = ? ######
non_zero_S <- 0
for(i in 1:1222){
  for(j in 1:1222){
    if(i<j){
      if(Res[[4]][i,j] > 0){
        non_zero_S = non_zero_S + 1
      }
    }
  }
}
non_zero_S  ## 915 ##
S_hat = Res[[4]];

### TRUE Positive ###
TP = 0;
for(i in 1:586){
  for(j in 587:1222){
    if(S_hat[i,j]!=0 && pol[i,j]!=0){
       TP = TP + 1; 
    }
  }
}
TP

### Obtain the Ad_hoc_list
Ad_hoc_list <- c();
for(i in 1:586){
  for(j in 587:1222){
    if(S_hat[i,j]!=0 && pol[i,j]!=0){
      Ad_hoc_list <- rbind(Ad_hoc_list,c(i,j))
    }
  }
}
sort(table(Ad_hoc_list))

library('igraph'); library("CINNA"); library('greed')
data(Blogs)
X=as.matrix(Blogs[[1]])
colnames(X)

### False Positive ###
FP = 0;
for(i in 1:586){
  for(j in 587:1222){
    if(S_hat[i,j]!=0 && pol[i,j]==0){
      FP = FP + 1; 
    }
  }
}
FP

### False Negative ###
FN  = 0;
for(i in 1:586){
  for(j in 587:1222){
    if(S_hat[i,j]==0 && pol[i,j]!=0){
      FN = FN + 1; 
    }
  }
}
FN

### True Negative ###
TN  = 0;
for(i in 1:586){
  for(j in 587:1222){
    if(S_hat[i,j]==0 && pol[i,j]==0){
      TN = TN + 1; 
    }
  }
}
TN

### Non-zero element of the upper triangular part of S_hat ###
non_zero = 0;
for(i in 1:586){
  for(j in 587:1222){
    if(pol[i,j]!=0){
      non_zero = non_zero + 1; 
    }
  }
}
