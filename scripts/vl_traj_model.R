suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(library(cmdstanr))
suppressMessages(require(HDInterval))
suppressMessages(source('scripts/utils.R'))


inv_logit = function(x){
	exp(x)/(1+exp(x))
}


get_multi_idx = function(x, idx){
	return(as.numeric(gsub('\\]', '', str_split(str_split(x, '\\[', simplify=TRUE)[,2], ",", simplify=TRUE)[,idx])))
}

get_idx = function(x, idx){
	return(as.numeric(gsub('\\]', '', str_split(x, '\\[', simplify=TRUE)[,2])))
}


generate_dm = function(d, model=c()){
  if (!all(is.na(model))){
    require(tidyverse)
    make_indicator_vars = function(d2, cov){
      return(d2 %>% 
            mutate(idx = seq(1,n()), x=1) %>%
            select(all_of(c('idx', cov, 'x'))) %>%
            pivot_wider(
              names_from=all_of(cov),
              values_from=x, 
              values_fill=0) %>%
            select(-idx) %>%
            rename_with(~paste0(cov, .), names(.)) %>%
            select(all_of(sort(names(.)))))
    }
    # create empty design matrix 
    dm = d[,0]
    # log the category of each design matrix
    beta_cats = c()
    n_beta_cats = 1
    for (cov_idx in 1:length(model)){
      cov = model[cov_idx]
      cov_dms = list()
      # split on "*" to account for interaction terms
      # assumes at most pairwise interaction, 
      # could probably adjust but no at present
      cov_split = str_split(cov, "\\*", simplify=TRUE)[1,]
      for (i_cov in cov_split){
        if (is.numeric(d[[i_cov]])){
          cov_dms[[length(cov_dms)+1]] = d[i_cov]
        }else{
          cov_dms[[length(cov_dms)+1]] = make_indicator_vars(d, i_cov)
        }
      }
      # in the case of no interaction term
      if (length(cov_dms) == 1){
        cov_dms[[2]] = as.matrix(rep(1,nrow(d)))
      }
      # iterate over levels of the *second* covariate and
      # generate a design matrix for the *first*
      # coeffs within the first covariate sum to zero given the second
      for (col in 1:ncol(cov_dms[[2]])){
        dm = bind_cols(dm,
          tibble(cov_dms[[1]] * t(cov_dms[[2]][,col])) %>%
            rename_with(~paste0(., "*", colnames(cov_dms[[2]])[col])))
        beta_cats = c(beta_cats, rep(n_beta_cats,ncol(cov_dms[[1]])))
        n_beta_cats = n_beta_cats + 1
      }
    }
    # remove trailing colon when no interaction term
    colnames(dm) = gsub("\\*$", "", colnames(dm))
    out=list()
    out$dm = dm
    out$beta_cats = beta_cats
  }else{
    out = list()
    out$dm = vector('numeric')
    out$beta_cats = vector('character')
  }
  return(out)
}


get_hpd = function(x){
	suppressMessages(require(HDInterval))
  return(
  	bind_cols(
  		enframe(hdi(x, credMass=0.95)) %>% pivot_wider(names_prefix="095"),
  		enframe(hdi(x, credMass=0.5)) %>% pivot_wider(names_prefix="050")))
}


#### MAIN ANALYSIS BEGINS HERE ####
labels = c('overall', 'nnrti', 'nrti')
dm_cols = c('round', 'age_cat', 'sex', 'comm_type')

d = read_tsv('data/rakai_drug_resistance_categorized_R15_R20.tsv', show_col_types=FALSE)  %>%
	group_by(round) %>%
	mutate(
		median_int_date = median(int_date, na.rm=TRUE)) %>%
	ungroup() %>%
	mutate(
		comm_type = case_when(
			comm_type == "Agrarian community" ~ "inland",
			comm_type == "Fish landing site" ~ "fishing",
			comm_type == "Trading community" ~ 'inland'),
		numeric_round = as.numeric(round),
		round = as.character(round))

# viremic PLHIV in round 16 thorugh 19
# with leading variables
# how many have na viremia status? 
ahead_d = d %>%
	filter(
		finalhiv == "P" & 
		!is.na(viremic) & 
		round >= 16) %>%
	arrange(round) %>%
	group_by(study_id) %>%
	mutate(
		lead_numeric_round = lead(numeric_round),
		lead_numeric_copies = lead(numeric_copies),
		lead_vl_copies_cat = lead(vl_copies_cat),
		lead_viremic = lead(viremic)) %>%
	ungroup() %>%
	group_by_at(dm_cols) %>%
	mutate(strata_size = sum(!is.na(lead_viremic))) %>%
	ungroup() %>%
	mutate(overall = 'susceptible') %>%
	filter(viremic == TRUE & round < 20)

# label results by the followup rounds
uniq_rounds = d %>% 
	filter(round >= 17 & round <= 20) %>%
	ungroup() %>%
	select(round, median_int_date) %>% 
	unique() %>% 
	arrange(round) %>% 
	mutate(
		round_idx = seq(1,n()),
		median_int_date = as.Date(median_int_date))

#### PROBABILITY OF SUPPRESSION AMONG VIREMIC PLHIV BY RESISTANCE STATUS ####
m = cmdstanr::cmdstan_model(paste0('stan/vl_pairs_resistance.stan'))
for (i in c(0, 1)){
	for (drug_class in c('overall', 'nnrti', 'nrti')){
		# get size of strata defined by dm_cols
		# then filter to those that participated (as identified by having VL) in the following round
		# and filter for those with NNRTI resistance data in this round
		# then add inidividual identifier
		use_d = ahead_d %>% 
			filter(
				!is.na(get(drug_class))) %>%
			group_by(study_id) %>%
			mutate(
				idx = cur_group_id(),
				lead_suppressed = case_when(
					is.na(lead_viremic) ~ -1,
					lead_numeric_round != numeric_round + 1 ~ -1,
					(lead_viremic == FALSE) & (lead_numeric_round == numeric_round + 1) ~ 1, 
					(lead_viremic == TRUE) & (lead_numeric_round == numeric_round + 1) ~ 0)) %>%
			group_by_at(dm_cols) %>%
			mutate(
				obs_strata_size_susc = sum(!!as.name(drug_class) == 'susceptible'),
				obs_strata_size_res = sum(!!as.name(drug_class) == 'intermediate/high')) %>%
			ungroup()
		dm = generate_dm(use_d, dm_cols)
		# using group by/slice so order is consistent when we get weights
		uniq_use_d = use_d %>%
				group_by_at(dm_cols) %>% 
				filter(row_number() == 1) %>% 
				select(all_of(c(dm_cols, 'strata_size'))) %>%
				ungroup()
		uniq_dm = generate_dm(uniq_use_d, dm_cols)
		stopifnot(all(colnames(dm$dm) == colnames(uniq_dm$dm)))
		stan_data = list(
					# input data characteristics
					N_ind = max(use_d$idx),
					N_pair = nrow(use_d),
					ind_idx = use_d$idx,
					# design matrix
					N_round = sum(grepl('round', colnames(dm$dm))),
					N_age = sum(grepl('age_cat', colnames(dm$dm))),
					N_sex = sum(grepl('sex', colnames(dm$dm))),
					N_comm = sum(grepl('comm_type', colnames(dm$dm))),
					X = as.matrix(dm$dm),
					# outcome data
					y = use_d$lead_suppressed,
					r = 1*(use_d[[drug_class]] == 'intermediate/high'),
					# run parameters
					sample_from_posterior = i,
					# supports generated quantities
					N_uniq = nrow(uniq_dm$dm),
					X_uniq = uniq_dm$dm,
					N_per_uniq = uniq_use_d$strata_size,
					N_per_uniq_obs = c(
						(use_d %>%
							group_by_at(dm_cols) %>% 
							filter(row_number() == 1))$obs_strata_size_susc,
						(use_d %>%
							group_by_at(dm_cols) %>% 
							filter(row_number() == 1))$obs_strata_size_res))
		fit = m$sample(
			    data = stan_data,
			    seed = 42, chains = 4, parallel_chains = 4,
			    iter_warmup = 2000, iter_sampling = 2000, refresh = 500)
		fit_draws = as_tibble(fit$draws(inc_warmup = FALSE, format = "draws_df")) %>%
			mutate(iter = seq(1,n()))
		fit_sum = as_tibble(bind_rows(
			fit_draws[grepl('alpha_0', colnames(fit_draws)) | 
					grepl('alpha_sd', colnames(fit_draws)) |
					grepl('betas\\[', colnames(fit_draws)) |
					grepl('mu_pred', colnames(fit_draws)) |
					grepl('rr', colnames(fit_draws))] %>%
				mutate(cit = seq(1,n())) %>%
				pivot_longer(-cit) %>%
			    group_by(name) %>% 
			    group_map(
					~cbind(get_hpd(.x$value) %>%
				    	mutate(median=median(.x$value)), .y)))) %>%
			mutate(
				drug_class=drug_class,
				posterior = as.logical(i),
				round_idx = case_when(
					# if item is two dimensional
					grepl('round_mu_pred', name) |
						grepl('between_round_rr', name) ~ get_multi_idx(name, 1),
					# if item is one dimensional
					grepl('within_round_rr', name) ~ get_idx(name)),
				status = case_when(
					# if item is two dimensional
					grepl('round_mu_pred', name) |
						grepl('between_round_rr', name) |
						grepl('strata_mu_pred', name) ~ c('susceptible', 'intermediate/high')[get_multi_idx(name, 2)],
					# if item is one dimensional
					grepl('overall_mu_pred', name) ~ c('susceptible', 'intermediate/high')[get_idx(name)])) %>%
			left_join(uniq_rounds, by='round_idx')
		# get and summarize parameters from draws that we want
		# save all
		saveRDS(fit, file = paste0('models/', drug_class, '_', c('prior', 'posterior')[i+1], '_prob_suppression.Rds'))
		write_tsv(fit_sum, paste0('models/', drug_class, '_', c('prior', 'posterior')[i+1], '_prob_suppression.tsv'))
	}
}

