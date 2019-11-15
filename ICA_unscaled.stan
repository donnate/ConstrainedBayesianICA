data{
	int N;  /// number of data points
	int T;  /// dimension (time)
	int<lower=1> K;  /// latent dimension
	matrix[T, N] X;  //data
	matrix[N, N] Structural_precision;  // regularization corresponding to structural information
	
	/// hyper parameters for the model

}

parameters{
   matrix[K, N] A_tilde;    // loadings
   vector[N] mu_A0; //means of the signals for A (will have to make them positive)
   
   positive_ordered[K] Lambda;
   matrix<lower=0, upper=1>[K, N] D;
   matrix[T, K] S;   // matrix sources latent signal
   //real<lower=0, upper = 1> pi_Lambda; //for variance of spike and slab
   vector<lower=0>[N] sigma_noise_inv;
}

transformed parameters{
	vector<lower=0>[N] sigma_noise;
	vector<lower=0>[N] mu_A; //means of the signals for A
	matrix<lower=0>[K, N] A;    // loadings
	
	sigma_noise= inv(sigma_noise_inv);
	mu_A = fabs(mu_A0);
	A = fabs(A_tilde);
}

model{
	for (j in 1:N){
			mu_A0[j]  ~ normal(0, 2);
	}
	
	for (k in 1:K){    //	     
		A_tilde[k,] ~ multi_normal_prec(mu_A, Structural_precision);
		Lambda[k]   ~ beta(5.0/K, 5.0/K);
		//target += log_mix(0.5, gamma_lpdf(Lambda[k]  | 1, 1), gamma_lpdf(Lambda[k]  | 10, 200));
		//Lambda[k] ~ gamma(K, K);
		for (j in 1:N){
			D[k, j]  ~  beta(0.5, 5);
			//D[k, j]  ~  beta(10.0/N, 1 - 10.0/N);
			//target += log_mix(0.2, gamma_lpdf(D[k, j]  | 2, 1), gamma_lpdf(D[k, j]  | 10, 100));
	    }
	}

	//// Generate the sources
	for (i in 1:K){
		for (t in 1:T){
			S[t,i] ~ double_exponential(0, 1/sqrt(2*T));
		}
	}
	
	//pi_Lambda ~ beta(5.0/K, 1- 5.0/K);
	for (i in 1:N){
             sigma_noise_inv[i] ~ exponential(1.0/10);
             to_vector(X[,i]) ~ normal(to_vector(S*diag_matrix(Lambda)* ((A .* D)[,i])), sigma_noise[i]);
	}
       
        //to_vector(X) ~ normal(to_vector(S*diag_matrix(Lambda)* (A .* D)), sigma_noise);
}

generated quantities{
	matrix[T, N] Y;
	matrix[K, N] Loadings;
	Y  = S * diag_matrix(Lambda) *  (A .* D);
	Loadings  = diag_matrix(Lambda) *  (A .* D);
}

