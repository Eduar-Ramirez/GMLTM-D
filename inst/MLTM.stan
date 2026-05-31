data {
  int N_subj;
  int N_item;
  array[N_subj * N_item] int ID;
  array[N_subj * N_item] int item;
  array[N_subj * N_item] int ones;
  int M;
  int K;
  int n_eta;
  array[n_eta, 2] int indexes_eta;
  matrix[N_item, K] Q;
  matrix[N_item, M] C;
  array[N_subj * N_item] int y;
  real prior_theta_mu;
  real<lower=0> prior_theta_sigma;
  real prior_eta_mu;
  real<lower=0> prior_eta_sigma;
  real prior_alpha_mu;
  real<lower=0> prior_alpha_sigma;
}
transformed data {
  int N = N_subj * N_item;
  matrix[N, M] Cx = C[item, ];
}
parameters {
  matrix<lower=0>[1, M] alpha;
  vector[n_eta] eta;
  matrix[N_subj, M] theta;
}
model {
  vector[N] p;
  matrix[N, M] mu;
  matrix[N_item, M] beta;
  matrix[K, M] eta_matrix;

  eta_matrix = rep_matrix(0, K, M);
  for(i in 1:n_eta) {
    eta_matrix[indexes_eta[i, 1], indexes_eta[i, 2]] = eta[i];
  }
  beta = Q * eta_matrix;
  mu = inv_logit(alpha[ones, ] .* (theta[ID, ] - beta[item, ])) .^ Cx;
  for(i in 1:N) {
    p[i] = prod(mu[i, ]);
  }
  target += normal_lpdf(to_vector(theta) | prior_theta_mu, prior_theta_sigma);
  target += normal_lpdf(to_vector(alpha) | prior_alpha_mu, prior_alpha_sigma) -
            normal_lccdf(0 | prior_alpha_mu, prior_alpha_sigma);
  target += normal_lpdf(eta | prior_eta_mu, prior_eta_sigma);
  target += bernoulli_lpmf(y | p);
}
generated quantities {
  vector[N] log_lik;
  {
    vector[N] p;
    matrix[N, M] mu;
    matrix[N_item, M] beta;
    matrix[K, M] eta_matrix;
    eta_matrix = rep_matrix(0, K, M);
    for(i in 1:n_eta) {
      eta_matrix[indexes_eta[i, 1], indexes_eta[i, 2]] = eta[i];
    }
    beta = Q * eta_matrix;
    mu = inv_logit(alpha[ones, ] .* (theta[ID, ] - beta[item, ])) .^ Cx;
    for(i in 1:N) {
      p[i] = prod(mu[i, ]);
      log_lik[i] = bernoulli_lpmf(y[i] | p[i]);
    }
  }
}
