GD <- function(para_v){
  l <- 1
  N <- para_v[1]
  lambda <- para_v[2]
  
  X_a <- para_v[3]
  Z_a <- para_v[4]
  U_a <- para_v[5]
  
  X_M <- matrix(para_v[6:((N^2)+5)],N,N)
  Z_M <- matrix(para_v[((N^2)+6):((2*(N^2)+5))],N,N)
  U_M <- matrix(para_v[((2*(N^2)+6)):((3*N^2)+5)],N,N)
  
  X <- matrix(para_v[((3*(N^2)+6)):((4*N^2)+5)],N,N)
  
  new_sol <- rep(0,(N^2+1))
  old_sol <- rep(0,(N^2+1))
  Hess_M = matrix(0, N, N)
  
  grad <- rep(0,(N^2+1))
  g <- matrix(0,N,N)
  K = -sum(X[upper.tri(X,diag=FALSE)])
  P = Z_a-U_a; Q = Z_M-U_M;
  
  while(TRUE){
    A = exp(X_a + X_M)/(1+exp(X_a + X_M))
    
    # Compute the alpha : gradient/hessian
    grad[1] = ( K + sum(A[upper.tri(A,diag=FALSE)]) + (1/lambda)*(X_a-P) )
    # grad[1] = ( K + (1/N)*sum(A) + (1/lambda)*(X_a-P) )
    
    # Compute the Gradient of M && Compute the Hessian of M
    g = (-1/2)*X + (1/lambda)*(X_M-Q)
    g[upper.tri(g,diag = FALSE)] = g[upper.tri(g,diag = FALSE)] + A[upper.tri(A,diag=FALSE)]
    # g = (-1/N)*X + (1/lambda)*(X_M-Q) + (1/N)*A
    grad[2:(N^2+1)] = as.vector(g)
    
    # Update the solution & Calculate the difference between objective function
    old_sol = c(X_a, X_M)
    new_sol = old_sol - 0.05*grad
    CC<-max(abs(grad))
    print(CC)
    # Convergence Criteria
    if(CC < 10^-7){
      X_a <- new_sol[1]
      X_M <- matrix(new_sol[2:(N^2+1)],N,N)
      break
    }
    # Store the variable
    X_a <- new_sol[1]
    X_M <- matrix(new_sol[2:(N^2+1)],N,N)
    l = l + 1
  }
  return(new_sol) 
}