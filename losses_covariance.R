##### LOSSES ##########
source('/scratch/users/cdonnat/Bayesian_ICA/auxiliary_functions.R')

operator_norm <- function(A, B){
  l = sort(abs(eigen(A-B, only.values = TRUE)$values))
  return(max(l))
}

norm_eigenbulk <- function(A, B){
  l_A = sort(eigen(A, only.values = TRUE)$values)
  l_B = sort(eigen(B, only.values = TRUE)$values)
  return(sum((l_A-l_B)^2))
}

norm_eigen_KS <- function(A, B){
  l_A = sort(eigen(A, only.values = TRUE)$values)
  l_B = sort(eigen(B, only.values = TRUE)$values)
  return(max(abs((l_A-l_B))))
}


loss_em <- function(Sigma_m1, Sigma_hat_m1, S, T){
  return(sum(diag((Sigma_m1- Sigma_hat_m1 )%*% (Sigma_m1- Sigma_hat_m1 ) %*%S))/(T * sum(diag(Sigma_m1))))
}


loss_n <- function(A,B){
  l_A = sort(abs(eigen(A-B, only.values = TRUE)$values))
  return(sum(sort(l_A)))
}

loss_stein <- function(A,B){
    #A_m1 = solve(A)
    #A_m1 = svd(A)
  dec = eigen(A, symmetric = TRUE)
  D <- diag(sapply(dec$values, FUN=function(x){ifelse(x>1e-7, 1/x, 0)}))
  A_m1 = dec$vectors %*% D %*% t(dec$vectors)
  N = nrow(A)
  return(sum(diag(A_m1%*%B - eye(N))) - 0.5* log( det(B)/max(10^(-7),det(A))))
  
}



loss_frechet <- function(A,B){
  svd_A = svd(A)
  sqrt_A = svd_A$u %*% diag(sqrt(svd_A$d))%*% t(svd_A$u)
  svd_B = svd(B)
  sqrt_B = svd_B$u %*% diag(sqrt(svd_B$d))%*% t(svd_B$u)
  return(sum(diag(A + B -2 * sqrt_A%*%sqrt_B)))
  
}

loss_frob <- function(A,B){
  return(sum((A-B)^2))
}


loss_fpr <- function(A,B){
  B_tilde = (A > 0 )
  return(sum((A <=0 ) * (B>0))/sum(B>0))
}

loss_tpr <- function(A,B){
  return(sum((A > 0 ) * (B>0))/sum(A>0))
}



stein_loss <- function(Sigma_hat, Sigma_inv) {
    if (any(is.na(Sigma_hat))) return(NA)
    return(sum(diag(Sigma_hat%*%Sigma_inv)) - log(det(Sigma_hat%*%Sigma_inv)) - nrow(Sigma_inv))
}
