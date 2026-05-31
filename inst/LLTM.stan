data {
  int N_subj;
  int N_item;
  array[N_subj*N_item] int ID;
  array[N_subj*N_item] int item;
  int K;
  matrix[N_item, K] Q;
  array[N_subj*N_item] int y;
  real prior_theta_mu;
  real<lower=0> prior_theta_sigma;
  real prior_eta_mu;
  real<lower=0> prior_eta_sigma;
}
transformed data {
  int N = N_subj*N_item;
}
parameters {
  vector[N_subj] theta;
  vector[K] eta;
}
transformed parameters {
  vector[N_item] beta;
  beta = Q * eta;
}
model {
  vector[N] p;
  p = inv_logit(theta[ID] - beta[item]);
  target += normal_lpdf(theta | prior_theta_mu, prior_theta_sigma);
  target += normal_lpdf(eta   | prior_eta_mu,   prior_eta_sigma);
  target += bernoulli_lpmf(y | p);
}
generated quantities {
  vector[N] log_lik;
  {
    vector[N] p;
    p = inv_logit(theta[ID] - beta[item]);
    for(i in 1:N) {
      log_lik[i] = bernoulli_lpmf(y[i] | p[i]);
    }
  }
}
