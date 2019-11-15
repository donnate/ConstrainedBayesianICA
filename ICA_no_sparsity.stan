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

}

transformed parameters{
        real<lower=0> sigma_noise;
        matrix<lower=0>[K, N] A;    // loadings

        sigma_noise= inv(sigma_noise_inv);
        A = fabs(A1);
}

model{

        for (k in 1:K){    //
            A1[k,] ~ multi_normal_prec(rep_vector(0, N), Structural_precision);
                target += log_mix(0.2, gamma_lpdf(Lambda[k]  | 10., 10.), gamma_lpdf(Lambda[k]| 0.01, 1.));
        }


        //// Generate the sources
        for (i in 1:K){
                for (t in 1:T){
                        S[t,i] ~ double_exponential(0, 1.0/sqrt(2));
                }
        }
    sigma_noise_inv ~ gamma(1., 1.);
        to_vector(X) ~ normal(to_vector(S * diag_matrix(Lambda)* A), sigma_noise);
}

generated quantities{
        matrix[T, N] Y;
        matrix[K, N] Loadings;
        Y  = S * diag_matrix(Lambda) *  A ;
        Loadings  = diag_matrix(Lambda) *  A;
}
