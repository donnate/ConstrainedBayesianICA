##### Pipeline for the Bayesian ICA model
source('/Users/cdonnat/Dropbox/Bayesian_ICA/auxiliary_functions.R')
source('/Users/cdonnat/Dropbox/Bayesian_ICA/losses_covariance.R')
source('/Users/cdonnat/Dropbox/Bayesian_ICA/graph_generation.R')


	
library(rmutil)
require(rstan)
library(igraph)
library(MASS)

#Sys.setenv(USE_CXX14 = 1)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#args = commandArgs(trailingOnly=TRUE)  ### Pass the seed + name of saved file where we want to save the results.
seed_i = 1#args[1]
set.seed(seed_i)



#### Generate simple ICA models
N <- 80
M <- 20
Time <- 300
K_GT <- 5
noise_level = 0.1
n0 = rep(15, K_GT)
G = 1
B = 1
K =20






### A. We start with a few networks


Graphs <- list()
Cov <- list()
VecCov <- matrix(0, G*B, N*(N-1)/2)

for (g in 1:G){
  Graphs[[(g-1)*B + 1]] = create_graph(N,K_GT, n0)
  for ( b in 1:B){
    if (b>1){
      Graphs[[(g-1)*B + b]] = scramble(Graphs[[(g-1)*B + 1]]$Adj,
                                       Graphs[[(g-1)*B + 1]]$Graph,
                                       noise_level,
                                       n0)
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
mu = rep(2, length(Graphs[[1]]$Sigma_i))
n1 <- c(n0, ncol(Graphs[[1]]$Sigma_i[[K_GT+1]]))


for (g in 1:G){
  for ( b in 1:B){
    Factors_blocks[[(g-1)*B + b]]<- matrix(0, length(Graphs[[(g-1)*B + b]]$Sigma_i)-1, N)
    a = 1
    list_links <- c()
    for (i in 1:(length(Graphs[[(g-1)*B + b]]$Sigma_i)-1)){
      if (i==1){
        indices = a:(a+n1[i]-1)
      }
      else{
        indices = c(list_links, a:(a+n1[i]-length(list_links )-1))
      }
      print(indices)
      Factors_blocks[[(g-1)*B + b]][i, indices] = mvrnorm(1, rep(mu[i], n1[i]), 
                                                          Graphs[[(g-1)*B + b]]$Sigma_i[[i]])
      a = a+ n1[i]-length(list_links)
      list_links <- c(list_links,a-1)
    }
    
    S[[(g-1)*B + b]] = matrix(rlaplace(Time * (K_GT), s = 1/sqrt(2*Time)), Time, K_GT)
    Y[[(g-1)*B + b]] = S[[(g-1)*B + b]] %*% Factors_blocks[[(g-1)*B + b]] + matrix(rnorm(Time * (N), mean = 0, sd=1/sqrt(Time)), 
                                                                                   Time, N)
  }
}



########## FIT THE MODEL #################
m <- stan_model(file="/Users/cdonnat/Dropbox/Bayesian_ICA/ICA_sunday.stan")
Y.vb <- list()
A1.vb <- list()
A_tilde1.vb <- list()
S1.vb <- list()
mu_a0.vb <- list()
D1.vb <- list()
sigma_noise.vb <- list()
sigma_noise_inv.vb <- list()
Lambda.vb <- list()

A1.stan <- list()
A_tilde1.stan <- list()
S1.stan <- list()
mu_a0.stan <- list()
Lambda.stan <- list()
Loadings.stan <- list()
D1.stan <- list()
sigma_noise.stan <- list()



for (g in 1:G){
  for ( b in 1:B){
  
  
    Structural_prior = Graphs[[(g-1)*B + b]]$Adj
    Structural_prior[Structural_prior!=0]=1
    alpha = 0.1
    Lap = diag(apply(Structural_prior, 1, sum)) - Structural_prior + alpha * diag(N)
    
    data2 <- list(N = N, T = Time, K= 20, X = Y[[(g-1)*B + b]], 
                  Structural_precision= Lap)
    stan.fit.vb <- vb(m, data = data2, algorithm = "meanfield")
    
    
    
    A_tilde1.vb[[(g-1)*B + b]] <- apply(rstan::extract(stan.fit.vb,"A_tilde")[[1]], c(2,3), mean)
    S1.vb[[(g-1)*B + b]] <- apply(rstan::extract(stan.fit.vb,"S")[[1]], c(2,3), mean)
    mu_a0.vb[[(g-1)*B + b]] <- apply(rstan::extract(stan.fit.vb,"mu_A0")[[1]], 2 , mean)
    Lambda.vb[[(g-1)*B + b]] <- apply(rstan::extract(stan.fit.vb,"Lambda")[[1]], 2, mean)
    D1.vb[[(g-1)*B + b]] <- apply(rstan::extract(stan.fit.vb,"D")[[1]],  c(2,3), mean)
    sigma_noise.vb[[(g-1)*B + b]] = mean(rstan::extract(stan.fit.vb,"sigma_noise")[[1]])
    sigma_noise_inv.vb[[(g-1)*B + b]] = mean(rstan::extract(stan.fit.vb,"sigma_noise_inv")[[1]])
    
    initf1 <- function() {
           list(A_tilde=A_tilde1.vb[[(g-1)*B + b]], 
           mu_A0 = mu_a0.vb[[(g-1)*B + b]],
           D=D1.vb[[(g-1)*B + b]], S=S1.vb[[(g-1)*B + b]],
           Lambda = Lambda.vb[[(g-1)*B + b]], 
           sigma_noise_inv=sigma_noise_inv.vb[[(g-1)*B + b]])
    } 
    
    
    stan.fit = sampling(m, data = data2, chains=4,  warmup = 500, init = initf1, iter = 1000,
                    refresh = 100)
    
    Loadings.stan[[(g-1)*B + b]] <- apply(rstan::extract(stan.fit,"Loadings")[[1]], c(2,3), mean)                
    A_tilde1.stan[[(g-1)*B + b]] <- apply(rstan::extract(stan.fit,"A_tilde")[[1]], c(2,3), mean)
    S1.stan[[(g-1)*B + b]] <- apply(rstan::extract(stan.fit,"S")[[1]], c(2,3), mean)
    mu_a0.stan[[(g-1)*B + b]] <- apply(rstan::extract(stan.fit,"mu_A")[[1]], 2 , mean)
    Lambda.stan[[(g-1)*B + b]] <- apply(rstan::extract(stan.fit,"Lambda")[[1]], 2, mean)
    D1.stan[[(g-1)*B + b]] <- apply(rstan::extract(stan.fit,"D")[[1]],  c(2,3), mean)
    sigma_noise.stan[[(g-1)*B + b]] = mean(rstan::extract(stan.fit,"sigma_noise")[[1]])                
    
  }
}



########## SAVE EVERYTHING #################
sv_file = paste('/Users/cdonnat/Dropbox/Bayesian_ICA/results_synthetic_experiments_ICA/', save_folder,'_image_session.RData', sep='')
save.image(sv_file) 


