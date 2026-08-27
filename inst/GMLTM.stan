data {
  int N_subj;
  int N_item;
  array[N_subj * N_item] int ID;
  array[N_subj * N_item] int item;
  int M;
  int K;
  int n_eta;
  int n_alpha;
  array[n_eta, 2] int indexes_eta;
  array[n_alpha, 3] int indexes_alpha;
  matrix[N_item, K] Q;
  matrix[N_item, M] C;
  array[N_subj * N_item] int y;
  int cf;
  real prior_theta_mu;
  real<lower=0> prior_theta_sigma;
  real prior_eta_mu;
  real<lower=0> prior_eta_sigma;
  real prior_alpha_mu;
  real<lower=0> prior_alpha_sigma;
  real<lower=0> prior_c_shape1;
  real<lower=0> prior_c_shape2;
  int<lower=0, upper=1> alpha_lognormal; // 0 = semi-Normal (truncated at 0, default), 1 = Log-Normal
}
transformed data {
  int max_n_alpha = max(indexes_alpha[, 3]);
  int N = N_subj * N_item;
  matrix[N, M] Cx = C[item, ];
}
parameters {
  matrix[N_subj, M] theta;
  vector[n_eta] eta;
  vector<lower=0>[max_n_alpha] alpha;
  vector<lower=0, upper=1>[N_item] c;
}
model {
  matrix[K, M] eta_matrix;
  matrix[N_item, M] alpha_matrix;
  matrix[N_item, M] beta_matrix;
  matrix[N, M] mu;
  vector[N] p;

  eta_matrix = rep_matrix(0, K, M);
  alpha_matrix = rep_matrix(0, N_item, M);
  for(i in 1:n_eta) {
    eta_matrix[indexes_eta[i, 1], indexes_eta[i, 2]] = eta[i];
  }
  for(i in 1:n_alpha) {
    alpha_matrix[indexes_alpha[i, 1], indexes_alpha[i, 2]] = alpha[indexes_alpha[i, 3]];
  }
  beta_matrix = Q * eta_matrix;
  mu = inv_logit(alpha_matrix[item, ] .* (theta[ID, ] - beta_matrix[item, ])) .^ Cx;
  for(i in 1:N) {
    p[i] = c[item[i]] + (1-c[item[i]]) * prod(mu[i, ]);
  }
  target += normal_lpdf(to_vector(theta) | prior_theta_mu, prior_theta_sigma);
  target += normal_lpdf(eta | prior_eta_mu, prior_eta_sigma);
  if (alpha_lognormal == 1) {
    target += lognormal_lpdf(alpha | prior_alpha_mu, prior_alpha_sigma);
  } else {
    target += normal_lpdf(alpha | prior_alpha_mu, prior_alpha_sigma) -
              normal_lccdf(0 | prior_alpha_mu, prior_alpha_sigma);
  }
  target += beta_lpdf(c | prior_c_shape1, prior_c_shape2);
  target += bernoulli_lpmf(y | p);
}
generated quantities {
  vector[N] log_lik;
  {
    matrix[K, M] eta_matrix;
    matrix[N_item, M] alpha_matrix;
    matrix[N_item, M] beta_matrix;
    matrix[N, M] mu;
    vector[N] p;
    eta_matrix = rep_matrix(0, K, M);
    alpha_matrix = rep_matrix(0, N_item, M);
    for(i in 1:n_eta) {
      eta_matrix[indexes_eta[i, 1], indexes_eta[i, 2]] = eta[i];
    }
    for(i in 1:n_alpha) {
      alpha_matrix[indexes_alpha[i, 1], indexes_alpha[i, 2]] = alpha[indexes_alpha[i, 3]];
    }
    beta_matrix = Q * eta_matrix;
    mu = inv_logit(alpha_matrix[item, ] .* (theta[ID, ] - beta_matrix[item, ])) .^ Cx;
    for(i in 1:N) {
      p[i] = c[item[i]] + (1-c[item[i]]) * prod(mu[i, ]);
      log_lik[i] = bernoulli_lpmf(y[i] | p[i]);
    }
  }
}
