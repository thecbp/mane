 data {
  int<lower=0> J;              // number of subjects
  int<lower=0> K;              // number of treatments
  int<lower=0> N;              // sample size over all subjects
  matrix[N, K] X;              // predictor matrix of treatments
  real y[N];                   // outcomes
  int<lower=0, upper=J> id[N]; // patient identifier {1, ..., J}
}

parameters {
  corr_matrix[K] Omega;        // hyperprior, correlation matrix
  vector<lower=0>[K] tau;      // hyperprior, scales for correlation matrix
  vector[K] Gamma;             // hyperprior, population-level coefficient mean
  vector[K] Delta;             // for non-centered parameterization
  real<lower=0> y_sigma;       // within-subject variance
}

transformed parameters {
  matrix[J, K] Betas;          // prior, individual-level coefficients
  
  // Non-centered parameterization for better posterior exploration
  // Each person gets their own betas, sampled from the same population dist
  for (j in 1:J) {
    Betas[j,] = Gamma' + Delta' * quad_form_diag(Omega, tau);
  }
  
}

model {
  
  // Hyperpriors
  
    // Using a decomposition because it is easier to assign priors to
    // the scalar and correlation components and is recommended by Stan
    tau ~ normal(0, 5);
    Omega ~ lkj_corr(2);
    for (k in 1:K) {
      Gamma[k] ~ normal(0, 5);
    }
  
  // Priors
    y_sigma ~ normal(0, 5);
    
    // Centered parameterization for better posterior exploration
    // for (j in 1:J) {
    //   Betas[j,] ~ multi_normal(Gamma, quad_form_diag(Omega, tau)); 
    // }
    
  // Random scaling to perturb covariance matrix in non-centered approach
  // Only need a single Delta because all parameters come from the same 
  // population distribution
    for (k in 1:K) {
      Delta[k] ~ normal(0, 1);
    }
  
  // Likelihood 
    for (n in 1:N) {
      y[n] ~ normal(X[n,] * Betas[id[n],]', y_sigma); 
    }
       
}
