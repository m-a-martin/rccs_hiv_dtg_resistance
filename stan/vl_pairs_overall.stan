data {
  int<lower=1> N_ind;   // number of donor isnvs
  int<lower=1> N_pair; // number of transmission pairs
  array[N_pair] int ind_idx; // index of individual
  int N_round; // number of survey rounds
  int N_age; // number of age categories
  int N_sex; // number of sex categories
  int N_comm; // number of community categories
  matrix[N_pair, N_round + N_age + N_sex + N_comm] X; // design matrix
  int<lower = 1> N_uniq; // number of unique rows of design matrix, for generated quantities
  matrix[N_uniq, N_round + N_age + N_sex + N_comm] X_uniq; // unique rows of design matrix, for generated quantities
  vector[N_uniq] N_per_uniq; // total number of participants in each strata
  vector[N_uniq] N_per_uniq_obs; // number of observations in each strata
  array[N_pair] int y; // outcome variable
  int sample_from_posterior; // whether to sample from the posterior or the prior
}


parameters {
  real<lower = 0> alpha_sd;
  real alpha_0;
  vector[N_ind] alpha_i;
  sum_to_zero_vector[N_round] betas_round;
  sum_to_zero_vector[N_age] betas_age;
  sum_to_zero_vector[N_sex] betas_sex;
  sum_to_zero_vector[N_comm] betas_comm;
}

transformed parameters{
  vector[N_pair] logit_mu;
  vector[N_round + N_age + N_sex + N_comm] betas;
  betas[1:N_round] = betas_round;
  betas[N_round + 1: N_round + N_age] = betas_age;
  betas[N_round + N_age + 1: N_round + N_age + N_sex] = betas_sex;
  betas[N_round + N_age + N_sex + 1: N_round + N_age + N_sex + N_comm] = betas_comm;
  logit_mu = alpha_0 + alpha_sd * alpha_i[ind_idx] + X*betas;
}

model {
  alpha_0 ~ normal(0,1);
  alpha_sd ~ cauchy(0, 1);
  alpha_i ~ normal(0,1);
  //betas ~ normal(0,2);
  betas_round ~ normal(0,1);
  betas_age ~ normal(0,1);
  betas_sex ~ normal(0,1);
  betas_comm ~ normal(0,1);
  if (sample_from_posterior == 1){
      y ~ bernoulli(inv_logit(logit_mu));
  }
  
}

generated quantities{
  vector[N_uniq] strata_mu_pred;
  vector[N_round] round_mu_pred;
  vector[N_round] round_mu_pred_obs;
  vector[N_round] between_round_rr;
  vector[N_round] between_round_rr_obs;
  real overall_mu_pred;
  real overall_mu_pred_obs;
  strata_mu_pred = inv_logit(alpha_0 + X_uniq * betas);
  overall_mu_pred = sum( strata_mu_pred .* N_per_uniq ) /
    sum(N_per_uniq);
  overall_mu_pred_obs = sum( strata_mu_pred .* N_per_uniq_obs ) /
    sum(N_per_uniq_obs);
  // assumes round are the first columns in X_uniq
  for (i in 1:N_round){
    round_mu_pred[i] = sum( inv_logit(alpha_0 + X_uniq * betas) .* N_per_uniq_obs .* X_uniq[:,i] ) /
      sum( N_per_uniq_obs .* X_uniq[:,i] );
    round_mu_pred_obs[i] = sum( inv_logit(alpha_0 + X_uniq * betas) .* N_per_uniq .* X_uniq[:,i] ) /
      sum( N_per_uniq .* X_uniq[:,i] );
    between_round_rr[i] = round_mu_pred[i] / round_mu_pred[1];
    between_round_rr_obs[i] = round_mu_pred_obs[i] / round_mu_pred_obs[1];
  }
}

