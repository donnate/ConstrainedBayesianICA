eye <- function(n){
  return(diag(rep(1,n)))
}

nz_vals <- function(A,tol=1e-8) {A[upper.tri(A) & abs(A) > tol]}
nz_count <- function(A,tol=1e-8) {length(nz_vals(A,tol))}
lseq <- function(from, to, length.out) exp(seq(log(from), log(to), length.out=length.out))


projsplx = function(v, b=1){
  v[v<0] = 0
  u = sort(v, decreasing = TRUE)
  sv = cumsum(u)
  
  rho = which(u > (sv-b) / (1:length(v)))
  theta = max(0, (sv[rho] - b) / rho )
  w = v - theta
  w[w<0] = 0
  return(w)
}


laplacian <- function(A){
    return(diag(apply(A, 1, sum)) - A)
}

reg_laplacian <- function(A){
    N = nrow(A)
    diag(A) <- 0
    diagonalA <-apply(A, 1, FUN = function(x){ifelse(sum(x) > 10^(-5), 1/sqrt(sum(x)), 0)})
    diagonalA <- diag(diagonalA)
    return(eye(N) - diagonalA %*% A %*% diagonalA)
}


make_correlation <- function(A){
    diagonalA <- diag(sqrt(1/diag(A)))
    temp = -diagonalA %*% A %*% diagonalA
    diag(temp) <-  1
    return(temp)
}



covar2cor <- function(Sigma){
    d = diag(1/sqrt(diag(Sigma)))
    return(d%*%Sigma%*%d)
    
}


covar2invcor <- function(Sigma){
    d = diag(1/sqrt(diag(Sigma)))
    return(solve(d%*%Sigma%*%d))
    
}
