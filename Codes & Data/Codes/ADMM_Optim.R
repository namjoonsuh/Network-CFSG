ADMM <- function(X, gamma, delta){
  library(RSpectra)
  N <- nrow(X)
  source('GD.R')
  # Initialization
  # Initialization on X^m = (X_a, X_M, X_L, X_S)
  X_a <- runif(1,0,1); X_M <- matrix(0,N,N);
  X_L <- matrix(0,N,N); X_S <- matrix(0,N,N);
  
  # Initialization on Z^m = (Z_a, Z_M, Z_L, Z_S)
  Z_a <- runif(1,0,1); Z_M <- matrix(0,N,N);
  Z_L <- matrix(0,N,N); Z_S <- matrix(0,N,N);
  
  # Initialization on U^m = (U_a, U_M, U_L, U_S)
  U_a <- runif(1,0,1); U_M <- matrix(0,N,N);
  U_L <- matrix(0,N,N); U_S <- matrix(0,N,N);
  
  J = diag(1,N)-(1/N)*matrix(1,N,N)
  
  lambda <- 0.5; k <- 1; r <- 4;
  Obj <- rep(0,500);
  
  #Obj_Val<-function(para_v){
  #    X_a <- para_v[1]
  #    X_M <- matrix(para_v[2:(N^2+1)],N,N)
    
  #    A <- log(1+exp(X_a + X_M))
  #    (-X_a/N)*sum(X[upper.tri(X,diag=FALSE)]) + (-1/(2*N))*sum(X*X_M) + 
  #      (1/N)*sum(A[upper.tri(A,diag=FALSE)]) + 
  #      (1/(2*lambda))*(X_a-(Z_a-U_a))^2 + (1/(2*lambda))*sum((X_M-(Z_M-U_M))^2) 
  #}
  
  for(k in 1:500){
    ########## Step 1. <- Update x = ( a, M, L, S ) ##########
    # a and M update via Gradient Descent
    para_v <- c(N, lambda, X_a, Z_a, U_a, X_M, Z_M, U_M, X)
    res <- GD(para_v)
    X_a <- res[1]
    X_M <- matrix(res[2:(N^2+1)],N,N)
    
    # a and M update via L-BFGS method
    #para <- c(X_a, X_M)
    #res <- optim(para, Obj_Val, method = "CG", lower = -Inf, upper = Inf)
    #X_a <- res[[1]][1]
    #X_M <- matrix(res[[1]][2:(N^2+1)],N,N)
    
    # L Update : Eigen-thresholding
    ev<-eigs_sym(Z_L-U_L,r); VAL <- ev$values; EV <- ev$vectors;
    eigen_threshold <- VAL-delta*lambda
    eigen_threshold[eigen_threshold<0]<-0
    D <- as.matrix(diag(eigen_threshold[1:r]))
    # print(D); #print(EV)
    X_L <- J%*%(EV[,1:r]%*%D%*%t(EV[,1:r]))%*%J
    
    #ev<-eigen(Z_L-U_L); EV <- ev$vectors;
    #eigen_threshold <- ev$values-(delta*lambda)
    #eigen_threshold[eigen_threshold<0]<-0
    #D <- diag(eigen_threshold)
    #X_L <- J%*%(EV%*%D%*%t(EV))%*%J
    
    # S Update : Soft_thresholding : sign(x)*max(|x|-gamma*lambda,0)
    for(i in 1:N){
      for(j in 1:N){
        if(i!=j){
          X_S[i,j] = sign(Z_S[i,j]-U_S[i,j])*max( abs(Z_S[i,j]-U_S[i,j])-gamma*lambda,0 )
        }
        if(i==j)
          X_S[i,j] = Z_S[i,j]-U_S[i,j]
      }
    }
    
    ## Step 2. <- Update z = ( Z_a, Z_M, Z_L, Z_S ) ########## 
    a_Bar <- X_a + U_a
    M_Bar <- X_M + U_M
    L_Bar <- X_L + U_L
    S_Bar <- X_S + U_S
    
    Z_a <- a_Bar
    Z_M <- (1/3)*M_Bar + (1/3)*t(M_Bar) + (1/3)*L_Bar + (1/3)*S_Bar
    Z_L <- (1/6)*M_Bar + (1/6)*t(M_Bar) + (2/3)*L_Bar - (1/3)*S_Bar  
    Z_S <- (1/6)*M_Bar + (1/6)*t(M_Bar) - (1/3)*L_Bar + (2/3)*S_Bar
    
    ## Step 3. <- Update u = ( U_M, U_L, U_S, U_a ) ##########
    U_a <- U_a + (X_a - Z_a)
    U_M <- U_M + (X_M - Z_M)
    U_L <- U_L + (X_L - Z_L)
    U_S <- U_S + (X_S - Z_S)
  
    # res <- (X_a-Z_a)^2 + sum((X_M-Z_M)^2) + sum((X_L-Z_L)^2) + sum((X_S-Z_S)^2)
    res <- sqrt(sum((X_M-X_L-X_S)^2))
    # print(res)
    # Convergence Check
    if(res < 1e-2){
      result <- list(X_a, X_M, X_L, X_S)
      break
    }
    k <- k + 1
  }
  
  result <- list(X_a, X_M, X_L, X_S)
  return(result)
} 