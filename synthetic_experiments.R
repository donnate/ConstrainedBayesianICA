##### Pipeline for the Bayesian ICA model
source('Bayesian_ICA/auxiliary_functions.R')
source('Bayesian_ICA/losses_covariance.R')
source('Bayesian_ICA/graph_generation.R')


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

Sys.setenv(USE_CXX14 = 1)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


args = commandArgs(trailingOnly=TRUE)  ### Pass the seed + name of saved file where we want to save the results.
seed_i = args[1]
set.seed(seed_i)



#### Generate simple ICA models
N <- 65
M <- 20
Time <- 300
K_GT <- 5
n0 = rep(15, K_GT)
G = 5
B = 1
K = 20

save_folder = args[2]
sig = as.numeric(args[3])
transform_prior = as.numeric(args[4])


### A. We start with a few networks


Graphs <- list()
Cov <- list()
VecCov <- matrix(0, G*B, N*(N-1)/2)
nodes2swap <- list()
for (g in 1:G){
    Graphs[[(g-1)*B + 1]] = create_connected_graph(N,K_GT, n0)
  for ( b in 1:B){
    if (b>1){
	    nodes2swap[[(g-1)*B + b]] = sample(1:length(V(Graphs[[(g-1)*B + 1]]$Graph)), 6, replace=FALSE)
	    Graphs[[(g-1)*B + b]] = scramble_connected(Graphs[[(g-1)*B + 1]], nodes2swap[[(g-1)*B + b]])
    }
    Cov[[(g-1)*B + b]] <- graph2cov(Graphs[[(g-1)*B + b]]$Adj)
    VecCov[(g-1)*B + b,] <- Cov[[(g-1)*B + b]]$Omega[upper.tri(Cov[[(g-1)*B + b]]$Omega, diag = FALSE)]
  }
}

Factors_blocks <- list()
Mixings <- list()
Factors <- list()
S <- list()
Y <- list()
mu = c(1.5,1.7,2,2.5,3)
n1 <- c(n0, ncol(Graphs[[1]]$Sigma_i[[K_GT]]))


for (g in 1:G){
  for ( b in 1:B){
    Factors_blocks[[(g-1)*B + b]]<- matrix(0, length(Graphs[[(g-1)*B + b]]$Sigma_i), N)
    a = 1
    list_links <- c()
    for (i in 1:(length(Graphs[[(g-1)*B + b]]$Sigma_i))){
      if (i==1){
        indices = a:(a+n1[i]-1)
      }
      else{
        indices = c(list_links, a:(a+n1[i]-length(list_links )-1))
      }
      print(indices)
      n_temp = nrow(Graphs[[(g-1)*B + b]]$Sigma_i[[i]])
      Factors_blocks[[(g-1)*B + b]][i, indices] = mvrnorm(1, rep(2, n_temp), 
                                                          Graphs[[(g-1)*B + b]]$Sigma_i[[i]])
      a = a+ n1[i]-length(list_links)
      list_links <- c(list_links,a-1)
    }
    
    S[[(g-1)*B + b]] = matrix(rlaplace(Time * (K_GT), s = 1/sqrt(2*Time)), Time, K_GT)
    Y[[(g-1)*B + b]] = S[[(g-1)*B + b]] %*% Factors_blocks[[(g-1)*B + b]] + matrix(rnorm(Time * (N), mean = 0, sd=sqrt(sig)), 
                                                                                   Time, N)
  }
}



########## FIT THE MODEL #################
m.full <- stan_model(file="ICA_final.stan")
m.ns <- stan_model(file="ICA_no_sparsity.stan")

Y.vb <- list()
A1.vb <- list()
A_tilde1.vb <- list()
S1.vb <- list()
mu_a0.vb <- list()
D1.vb <- list()
sigma_noise.vb <- list()
sigma_noise_inv.vb <- list()
Lambda.vb <- list()
Loadings.vb <- list()

A1.stan <- list()
A_tilde1.stan <- list()
S1.stan <- list()
mu_a0.stan <- list()
Lambda.stan <- list()
Loadings.stan <- list()
D1.stan <- list()
sigma_noise.stan <- list()

VanillaICA_A <-list()
VanillaICA_S <-list()
models <- list()
for (g in 1:G){
  for ( b in 1:B){
  
    ind = (g-1)*B + b
    Structural_prior = Graphs[[ind]]$Adj
    Structural_prior[Structural_prior!=0]=1
    alpha = 0.1
    # Lap = diag(apply(Structural_prior, 1, sum)) - Structural_prior + alpha * diag(N)
    diag_L = diag(1.0/sqrt(apply(Structural_prior, 1, sum)))
    Lap = (1 + alpha) * diag(N)  - diag_L %*% as.matrix(Structural_prior)  %*% diag_L
   if (transform_prior == 3){
      diag_L = diag(apply(Structural_prior, 1, sum))
      Lap = alpha * diag(N) + diag_L  - as.matrix(Structural_prior)
    }  

    data1 <- list(N = N, T = Time, K= 20, X = scale(Y[[ind]]), 
		  Structural_precision=diag(N))    
    data2 <- list(N = N, T = Time, K= 20, X = scale(Y[[ind]]), 
		  Structural_precision= Lap)

	  
    ### COMPARE WITH VANILLA ICA AS WELL
    A <- fastICA(scale(Y[[ind]]), 20, alg.typ = "parallel", alpha = 1,
                 method = "R", row.norm = FALSE, maxit = 200,
                 tol = 0.0001, verbose = TRUE)
   VanillaICA_A[[ind]] <- A$A
   VanillaICA_S[[ind]] <- A$S	  
   sv_file = paste('results_synthetic/FINAL_', save_folder, '_prior', toString(transform_prior), '.RData', sep='')
   save.image(sv_file) 
   old_A = abs(VanillaICA_A[[ind]])
   norm_A  = apply(old_A, 1, sum)
   perm = sort(norm_A,index.return=TRUE)$ix
   old_A = old_A[perm,]
   old_S = VanillaICA_S[[ind]][,perm]
   old_D = matrix(1,K, N)
   old_Lambda = norm_A[perm]/max(norm_A)
   old_Loadings= abs(VanillaICA_A[[ind]])
   old_sigma_noise_inv= 10
   initf1 <- function() {
        list(A=old_A, 
        D=old_D, S=old_S,
        Lambda = old_Lambda, 
        sigma_noise_inv=old_sigma_noise_inv)
    } 

    for (i in 1:4){
	    ind_res = (ind-1)* 4 + i 
	    if (i == 1){
		    data  = data2
		    m =  m.full 
		    
	    }
	    if (i == 2){
		    data  = data1
		    m =  m.full 
		    
	    }
	    if (i == 3){
		    data  = data2
		    m =  m.ns 
		    
	    }
	    if (i == 4){
		    data  = data1
		    m =  m.ns     
	    }
		
	    stan.fit = sampling(m, data = data, chains=1,  warmup = 800, init = initf1, iter = 2000,
			    refresh = 100)
	    #models[[ind_res]] <- stan.fit
            Loadings.stan[[ind_res]] <- apply(rstan::extract(stan.fit,"Loadings")[[1]], c(2,3), mean)                
	    A1.stan[[ind_res]] <- apply(rstan::extract(stan.fit,"A")[[1]], c(2,3), mean)
	    S1.stan[[ind_res]] <- apply(rstan::extract(stan.fit,"S")[[1]], c(2,3), mean)
	    Lambda.stan[[ind_res]] <- apply(rstan::extract(stan.fit,"Lambda")[[1]], 2, mean)
	    if (i < 3){
	        D1.stan[[ind_res]] <- apply(rstan::extract(stan.fit,"D")[[1]],  c(2,3), mean)
	    } else {
	        D1.stan[[ind_res]] <- matrix(0, K, N)
	    }
	    sigma_noise.stan[[ind_res]] = mean(rstan::extract(stan.fit,"sigma_noise")[[1]])
	    save.image(sv_file) 
    }
  }
}



