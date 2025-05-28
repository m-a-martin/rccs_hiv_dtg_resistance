data {
  int<lower=1> N_ind;   // number of viremic t0 individuals with resistance data
  int<lower=1> N_pair; // number of viremic t0 participant-visits with resistance data
  array[N_pair] int ind_idx; // index of individual
  int N_round; // number of survey rounds
  int N_age; // number of age categories
  int N_sex; // number of sex categories
  int N_comm; // number of community categories
  matrix[N_pair, N_round + N_age + N_sex + N_comm] X; // design matrix for all viremic t0 pairs with resistance data
  array[N_pair] int r; // t0 resistance status for each pair  (0 = suscpetible, 1 = resistance)
  array[N_pair] int y; // t1 viremia status (-1 = missing, 0 = viremic, 1 = suppressed)
  int<lower = 0, upper = 1> sample_from_posterior; // whether to sample from the posterior or the prior (0 = prior, 1 = posterior)
  int<lower = 1> N_uniq; // number of rows of unique design matrix, for generated quantities
  matrix[N_uniq, N_round + N_age + N_sex + N_comm] X_uniq; // unique rows of design matrix, for generated quantities
  vector[N_uniq] N_per_uniq; // total number of viremic participants in each strata at t0 (susceptible + resistant)
  vector[N_uniq*2] N_per_uniq_obs; // number of observed pairs in each strata, first as susceptible, then as resistant

}


transformed data{
  int N_pair_with_followup; // how many pairs have follow-up data
  N_pair_with_followup = 0;
  for (i in 1:N_pair){
    N_pair_with_followup += (y[i]==-1) ? 0 : 1; 
  }
  array[N_pair_with_followup] int N_pair_followup_idx; // index of pairs with follow-up data
  int j; 
  j = 1;
  for (i in 1:N_pair){
    if (y[i] >= 0){
      N_pair_followup_idx[j] = i;
      j += 1;
    }
  }
  matrix[N_pair, N_round*2] R; // matrix of interaction term between resistance and round
  matrix[N_uniq*2, N_round*2] R_uniq; // interaction term between resistance and round for unique strata,
                                      // once for susceptible, then for resistance
  // fill matrices with 0s
  R = rep_matrix(0, N_pair, N_round*2);
  R_uniq = rep_matrix(0, N_uniq*2, N_round*2);
  // assumes round is first columns in X
  // assumes each obs from exactly one round
  for (i in 1:N_round){
    R[:,(i-1)*2 + 1]  = X[:,i] .* (1 - to_vector(r)); 
    R[:,(i-1)*2 + 2]  = X[:,i] .* (to_vector(r));
    R_uniq[:N_uniq, (i-1)*2 + 1] = X_uniq[:,i] .* rep_vector(1, N_uniq);
    R_uniq[N_uniq+1:N_uniq*2,(i-1)*2 + 2] = X_uniq[:,i] .* rep_vector(1, N_uniq);
  }
  // whether to model resistance or not
  int<lower = 0, upper = 1>model_resistance; 
  model_resistance = ( sum(r) > 0 );
}

parameters {
  real<lower = 0> alpha_sd;
  real alpha_0;
  vector[N_ind] alpha_i;
  sum_to_zero_vector[N_round] betas_round;
  sum_to_zero_vector[N_age] betas_age;
  sum_to_zero_vector[N_sex] betas_sex;
  sum_to_zero_vector[N_comm] betas_comm;
  // hard coded number of rounds
  sum_to_zero_vector[2] R_betas1;
  sum_to_zero_vector[2] R_betas2;
  sum_to_zero_vector[2] R_betas3;
  sum_to_zero_vector[2] R_betas4;
  // for model of prob. resistance at t0
  // to do def prob_res_alpha_sd based on input criteria
  real prob_res_alpha_0;
  real prob_res_alpha_sd;
  vector[(model_resistance) ? N_ind: 0] prob_res_alpha_i;
  sum_to_zero_vector[N_round] prob_res_betas_round;
  sum_to_zero_vector[N_age] prob_res_betas_age;
  sum_to_zero_vector[N_sex] prob_res_betas_sex;
  sum_to_zero_vector[N_comm] prob_res_betas_comm;
}


transformed parameters{
  // for suppression model, conditional on follow-up
  vector[N_pair_with_followup] logit_mu; // logit probability of suppression among viremic t0 with resistance data
                                         // conditional on having follow-up
  vector[N_round + N_age + N_sex + N_comm] betas;
  vector[N_round*2] R_betas;
  betas[1:N_round] = betas_round;
  betas[N_round + 1: N_round + N_age] = betas_age;
  betas[N_round + N_age + 1: N_round + N_age + N_sex] = betas_sex;
  betas[N_round + N_age + N_sex + 1: N_round + N_age + N_sex + N_comm] = betas_comm;
  // hard coded number of rounds
  R_betas[1:2] = R_betas1;
  R_betas[3:4] = R_betas2;
  R_betas[5:6] = R_betas3;
  R_betas[7:8] = R_betas4;
  logit_mu = alpha_0 + 
    alpha_sd * alpha_i[ind_idx[N_pair_followup_idx]] + 
    X[N_pair_followup_idx,:]*betas + R[N_pair_followup_idx,:]*R_betas;
  // for t0 resistance model
  vector[N_pair] logit_pi; // logit probability of t0 resistance among viremic t0 with resistnace data
  vector[N_round + N_age + N_sex + N_comm] prob_res_betas;
  if (model_resistance){
    prob_res_betas[1:N_round] = prob_res_betas_round;
    prob_res_betas[N_round + 1: N_round + N_age] = prob_res_betas_age;
    prob_res_betas[N_round + N_age + 1: N_round + N_age + N_sex] = prob_res_betas_sex;
    prob_res_betas[N_round + N_age + N_sex + 1: N_round + N_age + N_sex + N_comm] = prob_res_betas_comm;
    logit_pi = prob_res_alpha_0 + prob_res_alpha_sd*prob_res_alpha_i[ind_idx] +  X*prob_res_betas;
  }
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
  // hard coded number of rounds
  R_betas1 ~ normal(0,1);
  R_betas2 ~ normal(0,1);
  R_betas3 ~ normal(0,1);
  R_betas4 ~ normal(0,1);
  // for model of resistnace among all viremic t0 participants
  prob_res_alpha_0 ~ normal(0,1);
  prob_res_alpha_sd ~ cauchy(0, 1);
  prob_res_alpha_i ~ normal(0,1);
  prob_res_betas_round ~ normal(0,1);
  prob_res_betas_age ~ normal(0,1);
  prob_res_betas_sex ~ normal(0,1);
  prob_res_betas_comm ~ normal(0,1);

  if (sample_from_posterior == 1){
    y[N_pair_followup_idx] ~ bernoulli(inv_logit(logit_mu));
    // only model resistance if there is at least one resistant sample
    // otherwise assume all resistant
    if (model_resistance){
      r ~ bernoulli(inv_logit(logit_pi));
    }
  }
}


generated quantities{
  vector[N_uniq] strata_pi_pred; // predicted resistance in each unique strata
  // for all, first column is susceptible, second is resistant
  matrix[N_uniq,2] N_per_uniq_pred; // number per uniq strata stratified by predicted resistance status
  matrix[N_uniq,2] strata_mu_pred; // predicted suppression in each epi strata defined by X_uniq
  vector[2] overall_mu_pred; // average predicted suppression weighted by total number of viremic t0 participants in each strata
  vector[2] overall_mu_pred_obs; // averaged predicted suppression weighted by number of observations in each strata
  matrix[N_round,2] round_mu_pred; //
  matrix[N_round,2] round_mu_pred_obs;
  matrix[N_round,2] between_round_rr;
  matrix[N_round,2] between_round_rr_obs;
  vector[N_round] within_round_rr;
  vector[N_round] within_round_rr_obs;
  // suppression by strata
  strata_mu_pred[:,1] = inv_logit(alpha_0 + X_uniq*betas + R_uniq[:N_uniq,:]*R_betas);
  if (model_resistance){
    strata_mu_pred[:,2] = inv_logit(alpha_0 + X_uniq*betas + R_uniq[N_uniq+1:2*N_uniq,:]*R_betas);
    // resistance by strata
    strata_pi_pred = inv_logit(prob_res_alpha_0 + X_uniq*prob_res_betas);
    N_per_uniq_pred[:,1] = N_per_uniq .* (1 - strata_pi_pred);
    N_per_uniq_pred[:,2] = N_per_uniq .* strata_pi_pred;
  }else{
    N_per_uniq_pred[:,1] = N_per_uniq;
    N_per_uniq_pred[:,2] = rep_vector(0, N_uniq);
  }


  // iterate over estimates for susceptible and resistant
  for (i in 1:(1 + model_resistance)){
    // overall suppression averaged over strata
    // first, weighted by number of participants in each
    overall_mu_pred[i] = sum( strata_mu_pred[:,i] .* N_per_uniq_pred[:,i] ) /
      sum(N_per_uniq_pred[:,i]);
    // weighted according to number of paired observations
    overall_mu_pred_obs[i] = sum( strata_mu_pred[:,i] .* N_per_uniq_obs[1 + (i - 1)*N_uniq:N_uniq*i]) /
      sum(N_per_uniq_obs[1 + (i - 1)*N_uniq:N_uniq*i]);
    // predicted suppression averaged over strata per round
    // assumes first columns of X are for round
    for (k in 1:N_round){
      // first, weighted by number of participants in each
      round_mu_pred[k,i] = sum( strata_mu_pred[:,i] .* N_per_uniq_pred[:,i] .* X_uniq[:,k] ) /
        sum( N_per_uniq_pred[:,i] .* X_uniq[:,k] );
      // weighted according to number of paired observations
      round_mu_pred_obs[k,i] = sum( strata_mu_pred[:,i] .*N_per_uniq_obs[1 + (i - 1)*N_uniq:N_uniq*i] .* X_uniq[:,k] ) /
        sum( N_per_uniq_obs[1 + (i - 1)*N_uniq:N_uniq*i] .* X_uniq[:,k] );
       // between round risk ratios
        between_round_rr[k,i] = round_mu_pred[k,i] / round_mu_pred[1,i];
        between_round_rr_obs[k,i] = round_mu_pred_obs[k,i] ./ round_mu_pred[1,i];
    }
  }
  if (model_resistance){
     // finally, within-round risk ratios
    within_round_rr = round_mu_pred[:,2] ./ round_mu_pred[:,1];
    within_round_rr_obs = round_mu_pred_obs[:,2] ./ round_mu_pred_obs[:,1];
  }
}
