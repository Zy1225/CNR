#' ---
#' title: Template code for running simulations involving spatially misaligned data; see Section 4 and Appendix C of the associated manuscript "Cokrig-and-Regress for Spatially Misaligned Data" for more details
#' abstract: Note the full simulations in the associated manuscript were run on multiple parallel sessions with different seeds. The below template code assumes the full simulations are run in a single session, and hence the results are subject to minor Monte Carlo error.
#' author: Originally written by ZYT
#' date: "Code started July 2023"
#' ---


##---------------------
#' # Load appropriate packages and files 
##---------------------
rm(list=ls())
library(metaSEM)
library(geoR)
library(sf)
library(Matrix)
library(spaMM)
library(splines)
here::i_am("Simulation/Template_Simulation_Script.R")
library(here)
source(here("Code","simulate_misaligned.R"))
source(here("Code","cnr.R"))
source(here("Code","bootstrap_cnr.R"))
source(here("Code","lnmr.R"))
source(here("simulation","simulation_result.R"))

#File containing location data for the response and covariates, together with the lower triangular part of the great-circle distance matrix constructed based on all locations
load(here("simulation","sim_loc.Rdata"))


##---------------------
#' ## Setting up
##---------------------
all_loc = rbind(sim_loc$x_loc,sim_loc$y_loc)
all_loc_sf = st_as_sf(all_loc, coords = c("lon", "lat"), crs = "WGS84", agr = "constant")
sim_dist_mat =  st_distance(all_loc_sf)
sim_dist_vec = sim_dist_mat[lower.tri(sim_dist_mat)] / 1000
sim_dist_vec = as.numeric(sim_dist_vec)
sim_dist_mat_y_x = as.matrix(st_distance(st_as_sf(sim_loc$y_loc, coords = c("lon", "lat"), crs = "WGS84", agr = "constant"),
                             st_as_sf(sim_loc$x_loc, coords = c("lon", "lat"), crs = "WGS84", agr = "constant"))
)


set.seed(201)
n_sim = 400
T_bootstrap = 250

##---------------------
#' ## List to store simulation results for each run
##---------------------
simulation_list = vector("list", n_sim)

##---------------------
#' ## Initialize true parameter values
##---------------------
K = 5
true_nu_vec_x = rep(0.5,K)
true_alpha_vec_x = rep(0.0015, K)
true_sigma2_vec_x = rep(1,K)
true_tau_vec_x = rep(0.15,K)
true_mu_vec_x = rep(0,K)

true_R <- outer(1:5, 1:5, FUN = "-") %>% 
  abs %>% 
  {0.5^{.}}
true_R[1,-1] = true_R[-1,1] = 0
true_R[2,-2] = true_R[-2,2] = 0

true_beta = c(2,1,0.5,1,0.5,1)
true_sigma2_rho = 0.2
true_nu_rho = 0.5
true_alpha_rho = 0.0015
true_tau_epsilon = 0.01

true_y_name = 'y'
true_x_name = paste0('X',1:K)
true_y_mean_model_formula = as.formula('~ X1 + X2 + X3 + X4 + X5')

for(s in 1:n_sim){
  ##---------------------
  #' ## Simulate spatially misaligned data
  ##---------------------
  sim_data = simulate_misaligned(nu_vec_x = true_nu_vec_x, alpha_vec_x = true_alpha_vec_x, sigma2_vec_x = true_sigma2_vec_x, tau_vec_x = true_tau_vec_x, mu_vec_x = true_mu_vec_x, R_mat = true_R, beta_vec = true_beta, sigma2_rho = true_sigma2_rho, nu_rho = true_nu_rho, alpha_rho = true_alpha_rho, tau_epsilon = true_tau_epsilon, y_loc = sim_loc$y_loc, x_loc =  sim_loc$x_loc, y_name = true_y_name, x_name = true_x_name, 
                                 y_mean_model_formula =  true_y_mean_model_formula, spline_basis_list = NULL, cnr_step3_model_fit = NULL, gc_dist = TRUE,
                                 bootstrap_use = FALSE, tchol_Sigma_mat = NULL, tchol_Sigma_err = NULL)
  
  ##---------------------
  #' ## Proposed CNR
  ##---------------------
  cnr_out = cnr(data_x = sim_data$data_x, 
                data_y = sim_data$data_y,  
                gc_dist = TRUE, 
                nu_vec_x = true_nu_vec_x, 
                dist_vec = sim_dist_vec, 
                y_name = true_y_name, 
                x_name = true_x_name, 
                loc_name = sim_data$loc_name,
                slmm_formula = as.formula(paste(true_y_name,' ~ 1  +',
                                                paste(true_x_name,collapse = '+'), ' + Matern(1|lon+lat)')), 
                nu_rho = true_nu_rho, 
                smooth_list = NULL, discretized_num = NULL, smooth_CI_level = NULL)
  
  ##---------------------
  #' ## Preliminary Bootstrap
  ##---------------------
  pre_bs = bootstrap_cnr(T_bootstrap = T_bootstrap, bootstrap_CI_level = 0.95, cnr_out = cnr_out, bootstrap_type = 'Prelim', pre_bs = NULL)
  
  ##---------------------
  #' ## Second Bootstrap
  ##---------------------
  second_bs = bootstrap_cnr(T_bootstrap = T_bootstrap, bootstrap_CI_level = 0.95, cnr_out = cnr_out, bootstrap_type = 'Second', pre_bs = pre_bs)
  
  ##---------------------
  #' ## Non-cross-correlated Bootstrap
  ##---------------------
  NCC_bs = bootstrap_cnr(T_bootstrap = T_bootstrap, bootstrap_CI_level = 0.95, cnr_out = cnr_out, bootstrap_type = 'NCC', pre_bs = pre_bs)
  
  ##---------------------
  #' ## 1-NMR
  ##---------------------
  lnmr_L1 = lnmr(l = 1, dist_mat_y_x = sim_dist_mat_y_x,
                 data_x = sim_data$data_x, 
                 data_y = sim_data$data_y, 
                 gc_dist = TRUE, 
                 dist_vec = sim_dist_vec,
                 y_name = true_y_name, 
                 x_name = true_x_name, 
                 loc_name = sim_data$loc_name,
                 slmm_formula = as.formula(paste(true_y_name,' ~ 1  +',
                                                 paste(true_x_name,collapse = '+'), ' + Matern(1|lon+lat)')), 
                 nu_rho = true_nu_rho, 
                 smooth_list = NULL, discretized_num = NULL, smooth_CI_level = NULL)
  
  ##---------------------
  #' ## 3-NMR
  ##---------------------
  lnmr_L3 = lnmr(l = 3, dist_mat_y_x = sim_dist_mat_y_x,
                 data_x = sim_data$data_x, 
                 data_y = sim_data$data_y,  
                 gc_dist = TRUE, 
                 dist_vec = sim_dist_vec, 
                 y_name = true_y_name, 
                 x_name = true_x_name, 
                 loc_name = sim_data$loc_name,
                 slmm_formula = as.formula(paste(true_y_name,' ~ 1  +',
                                                 paste(true_x_name,collapse = '+'), ' + Matern(1|lon+lat)')), 
                 nu_rho = true_nu_rho, 
                 smooth_list = NULL, discretized_num = NULL, smooth_CI_level = NULL)
  
  ##---------------------
  #' ## 5-NMR
  ##---------------------
  lnmr_L5 = lnmr(l = 5, dist_mat_y_x = sim_dist_mat_y_x,
                 data_x = sim_data$data_x, 
                 data_y = sim_data$data_y, 
                 gc_dist = TRUE, 
                 dist_vec = sim_dist_vec, 
                 y_name = true_y_name, 
                 x_name = true_x_name, 
                 loc_name = sim_data$loc_name,
                 slmm_formula = as.formula(paste(true_y_name,' ~ 1  +',
                                                 paste(true_x_name,collapse = '+'), ' + Matern(1|lon+lat)')), 
                 nu_rho = true_nu_rho, 
                 smooth_list = NULL, discretized_num = NULL, smooth_CI_level = NULL)
  
  simulation_list[[s]] = list(
    sim_data = sim_data, cnr_out = cnr_out, pre_bs = pre_bs, second_bs = second_bs, NCC_bs = NCC_bs, lnmr_L1 = lnmr_L1, lnmr_L3 = lnmr_L3, lnmr_L5 = lnmr_L5
  )
}


##---------------------
#' ## Process simulation results and assess performance
##---------------------
sim_out = simulation_result(simulation_list, true_beta = true_beta, true_sigma2_rho = true_sigma2_rho, true_alpha_rho = true_alpha_rho, true_tau_epsilon =  true_tau_epsilon)

sim_out$bias_beta_df
sim_out$rmse_beta_df
sim_out$ase_esd_ratio_df
sim_out$bias_theta_err_df
sim_out$rmse_theta_err_df
sim_out$beta_CI_coverage_df
sim_out$beta_CI_average_width_df
