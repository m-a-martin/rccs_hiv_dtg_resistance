#!/usr/bin/env bash

#### ----------------- ####
#### MAIN TEXT FIGURES ####
#### ----------------- ####
Rscript scripts/plot_hiv_prev.R
Rscript scripts/plot_resistance_among_viremic_tx.R
Rscript scripts/plot_resistance_among_viremic_pt.R
Rscript scripts/plot_prob_suppression.R
Rscript scripts/plot_inS153Y_dists.R

#### --------------------- ####
#### SUPPLEMENTARY FIGURES ####
#### --------------------- ####
Rscript scripts/plot_dtg_by_clinic.R
Rscript scripts/plot_inS153Y_freq.R
Rscript scripts/plot_inS153Y_co.R