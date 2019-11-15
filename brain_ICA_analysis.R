########################################################################################## 
############################## Load Libraries     ########################################
########################################################################################## 


is.installed <- function(mypkg){
    is.element(mypkg, installed.packages(lib.loc="~/R_libs")[,1])
  } 

if (!is.installed("rmutil")) {
	install.packages('rmutil', lib="~/R_libs")
	}
	
library(rmutil, lib.loc = "~/R_libs")
require(rstan)
library(igraph, lib.loc = "~/R_libs")
library(MASS)
library(fastICA, lib.loc="~/R_libs")


#Sys.setenv(USE_CXX14 = 1)
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())

args = commandArgs(trailingOnly=TRUE)  ### Pass the seed + name of saved file where we want to save the results.
sub = args[1]
session = args[2]
data_ts = args[3]
data_struct = args[4]
transform_ts = as.integer(args[5])
scale = as.integer(args[6])
scale = ifelse(scale==0, FALSE, TRUE)
save_dir = args[7]
prior_lowerdiag =  as.integer(args[8])
transform_prior =  as.integer(args[9])
readheader =  as.integer(args[10])
alpha = as.numeric(args[11])


########################################################################################## 
############################## Load the Data      ########################################
########################################################################################## 

if (readheader == 1){
  timeseries <- data.frame(t(read.table(data_ts, header = 1, row.names = 1,sep = ',')))
} else {
  timeseries <- data.frame(t(read.table(data_ts, header = FALSE,sep = ',')))
}

#print(timeseries)
timeseries[is.na(timeseries)] = rnorm(length(which(is.na(timeseries))))      #### deals with Missing data for certain voxels
if (transform_ts == 1){
   timeseries = log(-min(timeseries) + 1 + timeseries)
}
print(sum(is.na(timeseries)))

if (readheader == 1){
   struct <- data.frame(t(read.table(data_struct, header = 1, row.names = 1,sep = ',')))
} else{
   struct <- data.frame(t(read.table(data_struct, header = FALSE,sep = ' ')))
}

if (prior_lowerdiag ==1){
   D = diag(diag(as.matrix(struct)))
   diag(struct) <- 0 
   struct =  struct + t(struct) + D
}
if (transform_prior == 1){
  Structural_prior = data.frame(matrix(0, nrow(struct), ncol(struct)))
  Structural_prior[which(struct!=0, arr.ind = TRUE)] = 1    #### Converts our structural rmtrix prior into a binary adjacency matrix
} else{
   if (transform_prior == 2){
       Structural_prior = log(1+ struct)
   }else{
       Structural_prior <- struct
   }
}

N <- ncol(timeseries)
Time <- nrow(timeseries)
K <- 20

########################################################################################## 
############################## FIT BayesICA       ########################################
########################################################################################## 

inside= sqrt(apply(Structural_prior, 1, sum))
inv_inside  = sapply(inside, FUN=function(x){ ifelse(x<0.0001,0,1.0/x)}) 
diag_L = diag(inv_inside)
Lap = (1 + alpha) * diag(N)  - diag_L %*% as.matrix(Structural_prior)  %*% diag_L
if (transform_prior == 3){
   diag_L = diag(apply(Structural_prior, 1, sum))
   Lap = alpha * diag(N) + diag_L  - as.matrix(Structural_prior)
}

print(c(alpha, N, K, Time, min(Lap), min(timeseries)))
if (scale){
  m <- stan_model(file="ICA_final.stan")
  data2 <- list(N = N, T = Time, K= 20, X = scale(timeseries), Structural_precision= Lap)
} else{
  m <- stan_model(file="ICA_unscaled.stan")
  data2 <- list(N = N, T = Time, K= 20, X =  timeseries, Structural_precision= Lap)
}




### COMPARE WITH VANILLA ICA AS WELL
A <- fastICA(scale(timeseries), 20, alg.typ = "parallel", alpha = 1,
	 method = "R", row.norm = FALSE, maxit = 200,
	 tol = 0.0001, verbose = TRUE)
VanillaICA_A <- A$A
VanillaICA_S <- A$S	

old_A = abs(VanillaICA_A)
norm_A  = apply(old_A, 1, sum)
perm = sort(norm_A,index.return=TRUE)$ix
old_A = old_A[perm,]
old_S = VanillaICA_S[,perm]
old_D = matrix(1, K, N)
old_Lambda = norm_A[perm]/max(norm_A)
old_Loadings= abs(VanillaICA_A)
old_sigma_noise_inv= 10
initf1 <- function() {list(A=old_A, 
			D=old_D, S=old_S,
			Lambda = old_Lambda,Lambda1 = old_Lambda,
			sigma_noise_inv=old_sigma_noise_inv)
			} 


	  
  
########################################################################################## 
##### Sample using MCMC

stan.fit = sampling(m, data = data2, chains=1,  warmup = 1000, init = initf1, iter = 2000,
				refresh = 100)

if (scale){
   sv_file = paste(save_dir, '/final_D_rightsize_tighter_prior_scaled_element_wise_multiplication', toString(sub),'_session', toString(session), '_tts', toString(transform_ts), '_tp', toString(transform_prior), '_mcmc.RData', sep='')
} else{
   sv_file = paste(save_dir, '/final_D_rightsize_unscaled_element_wise_multiplication', toString(sub),'_session', toString(session), '_tts', toString(transform_ts), '_tp', toString(transform_prior), '_mcmc.RData', sep='')
}
save.image(sv_file)

Loadings.stan <- apply(rstan::extract(stan.fit,"Loadings")[[1]], c(2,3), mean)                
A1.stan <- apply(rstan::extract(stan.fit,"A")[[1]], c(2,3), mean)
S1.stan <- apply(rstan::extract(stan.fit,"S")[[1]], c(2,3), mean)
D1.stan <- apply(rstan::extract(stan.fit,"D")[[1]],  c(2,3), mean)
Lambda.stan <- apply(rstan::extract(stan.fit,"Lambda")[[1]], 2, mean)
sigma_noise.stan = mean(rstan::extract(stan.fit,"sigma_noise")[[1]])  


               
A1.stan.var <- apply(rstan::extract(stan.fit,"A")[[1]], c(2,3), var)
S1.stan.var <- apply(rstan::extract(stan.fit,"S")[[1]], c(2,3), var)
Lambda.stan.var <- apply(rstan::extract(stan.fit,"Lambda")[[1]], 2, var)
D1.stan.var <- apply(rstan::extract(stan.fit,"D")[[1]],  c(2,3), var)
sigma_noise.stan.var = var(rstan::extract(stan.fit,"sigma_noise")[[1]])  

########################################################################################## 
##############################   SAVE EVERYTHING ######################################### 
save.image(sv_file)
