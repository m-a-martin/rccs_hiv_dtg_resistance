#!/usr/bin/env bash

# format data
# not all input files shared publicly for privacy reasons
Rscript scripts/consolidate_rccs_dat.R
Rscript scripts/format_dr_dat.R

# sequence data analysis
# not all input files shared publicly for privacy reasons
bash scripts/process_seqs.sh

# DTG scale-up data
Rscript scripts/format_dtg_scaleup.R

# prevalence estimates
Rscript scripts/calc_prev.R
Rscript scripts/calc_mut_prev.R
Rscript scripts/calc_dtg_prev.R

# individual level trajectories
Rscript scripts/vl_traj_model.R

# tabulate data
bash scripts/tabulate_data.sh

# plot results
bash scripts/plot_results.sh