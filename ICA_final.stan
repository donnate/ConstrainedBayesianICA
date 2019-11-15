data {
	int N;  /// number of data points
	int T;  /// dimension (time)
	int<lower=1> K;  /// latent dimension
	matrix[T, N] X;  //data
	matrix[N, N] Structural_precision;  // regularization corresponding to structural information
	
	/// hyper parameters for the model

}

parameters{
   matrix[K, N] A1;    // loadings
   positive_ordered[K] Lambda;
   matrix[T, K] S;   // matrix sources latent signal
   real<lower=0> sigma_noise_inv;
   simplex[K] Dvec[N];
}

transformed parameters{
	real<lower=0> sigma_noise;
	matrix<lower=0>[K, N] A;    // loadings
	matrix<lower=0, upper=1>[K, N] D;
	sigma_noise= inv(sigma_noise_inv);
	for (i in 1:N){  
		D[,i] = sqrt(Dvec[i]);
	}
	A = fabs(A1);
}

model{
	for (k in 1:K){    //	 
	    A1[k,] ~ multi_normal_prec(rep_vector(0, N), Structural_precision);  
		target += log_mix(0.2, gamma_lpdf(Lambda[k]  | 10., 10.), gamma_lpdf(Lambda[k]| 0.01, 1.));
	}
	
	for (i in 1:N){
		Dvec[i] ~ dirichlet(rep_vector(1.0/(K), K));
	}

	//// Generate the sources
	for (i in 1:K){
		for (t in 1:T){
			S[t,i] ~ double_exponential(0, 1.0/sqrt(2));
		}
	}

	for (i in 1:N){
             sigma_noise_inv[i] ~ exponential(1.0/10);
             to_vector(X[,i]) ~ normal(to_vector(S*diag_matrix(Lambda)* ((A .* D)[,i])), sigma_noise[i]);
	}
       
}


generated quantities{
	matrix[T, N] Y;
	matrix[K, N] Loadings;
	Y  = S * diag_matrix(Lambda) * (A .*D);
	Loadings  = diag_matrix(Lambda) * (A .* D);
}