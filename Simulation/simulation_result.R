#' @title 
#' Process results for simulation study.
#' 
#' @description 
#' Process results for simulation study. See Section 4 and Appendix C of the associated manuscript "Cokrig-and-Regress for Spatially Misaligned Data" for more details.
#' 
#' @param simulation_list A list with each element representing the results of each simulation run. In particular, each element is a list with elements named 'sim_data' (returned by \code{simulate_misaligned} function), 'cnr_out' (returned by \code{cnr} function), 'pre_bs' (returned by \code{bootstrap_cnr} function with \code{bootstrap_type = 'Prelim'}), 'second_bs' (returned by \code{bootstrap_cnr} function with \code{bootstrap_type = 'Second'}), 'NCC_bs' (returned by \code{bootstrap_cnr} function with \code{bootstrap_type = 'NCC'}), 'lnmr_L1' (returned by \code{lnmr} function with \code{l=1}), 'lnmr_L3' (returned by \code{lnmr} function with \code{l=3}), and 'lnmr_L5' (returned by \code{lnmr} function with \code{l=5}).
#' @param true_beta A vector of true mean regression coefficients.
#' @param true_sigma2_rho True variance parameter of spatial random effect.
#' @param true_alpha_rho True Mat\'{e}rn range parameter for the spatial random effect.
#' @param true_tau_epsilon True residual variance parameter, or equivalently, true nugget parameter for the response.
#' 
#' @return A list with the elements being input arguments of the function call and the following additional elements:
#' \item{bias_beta_df:}{A data frame describing the empirical biases of different mean regression coefficient estimates, not including the intercept.}
#' \item{rmse_beta_df:}{A data frame describing the empirical root mean squared errors of different mean regression coefficient estimates, not including the intercept.}
#' \item{ase_esd_ratio_df:}{A data frame describing the ratio of average estimated standard error to the empirical standard deviation of mean regression coefficient estimates.}
#' \item{bias_theta_err_df:}{A data frame describing the empirical biases of different estimates for spatial random effect's variance parameter, spatial random effect's Mat\'{e}rn range parameter, and residual variance parameter.}
#' \item{rmse_theta_err_df:}{A data frame describing the empirical root mean squared errors of different estimates for spatial random effect's variance parameter, spatial random effect's Mat\'{e}rn range parameter, and residual variance parameter.}
#' \item{beta_CI_coverage_df:}{A data frame describing the empirical coverage probability of different 95% confidence intervals for the mean regression coefficients, not including the intercept.}
#' \item{beta_CI_average_width_df}:{A data frame describing the average width of different 95% confidence intervals for the mean regression coefficients, not including the intercept.}

simulation_result = function(simulation_list, true_beta, true_sigma2_rho, true_alpha_rho, true_tau_epsilon){
  
  #Extract beta point estimates
  hatbeta_cnr_mat = t(sapply(simulation_list, FUN = function(xlist){
    xlist$cnr_out$hatbeta
  }))
  
  hatbeta_bc_mat = t(sapply(simulation_list,function(xlist){
    2*xlist$cnr_out$hatbeta - xlist$pre_bs$bootstrap_hatbeta_mean
  }))
  
  hatbeta_pre_bs_mat = t(sapply(simulation_list, FUN = function(xlist){
    xlist$pre_bs$bootstrap_hatbeta_mean
  }))
  
  hatbeta_second_bs_mat = t(sapply(simulation_list, FUN = function(xlist){
    xlist$second_bs$bootstrap_hatbeta_mean
  }))
  
  hatbeta_NCC_bs_mat = t(sapply(simulation_list, FUN = function(xlist){
    xlist$NCC_bs$bootstrap_hatbeta_mean
  }))
  
  hatbeta_1nmr_mat = t(sapply(simulation_list, FUN = function(xlist){
    xlist$lnmr_L1$hatbeta
  }))
  
  hatbeta_3nmr_mat = t(sapply(simulation_list, FUN = function(xlist){
    xlist$lnmr_L3$hatbeta
  }))
  
  hatbeta_5nmr_mat = t(sapply(simulation_list, FUN = function(xlist){
    xlist$lnmr_L5$hatbeta
  }))
  
  #Bias of beta
  bias_beta_mat = rbind(apply(hatbeta_cnr_mat,2,mean) - true_beta,
                        apply(hatbeta_bc_mat,2,mean) - true_beta,
                        apply(hatbeta_1nmr_mat,2,mean) - true_beta,
                        apply(hatbeta_3nmr_mat,2,mean) - true_beta,
                        apply(hatbeta_5nmr_mat,2,mean) - true_beta)
  
  #RMSE of beta
  var_beta_mat = rbind(apply(hatbeta_cnr_mat,2,var),
                       apply(hatbeta_bc_mat,2,var),
                       apply(hatbeta_1nmr_mat,2,var),
                       apply(hatbeta_3nmr_mat,2,var),
                       apply(hatbeta_5nmr_mat,2,var))
  
  rmse_beta_mat = sqrt(bias_beta_mat^2 + var_beta_mat)
  
  #ASE/ESD of beta
  ase_esd_ratio_mat = rbind(
    #Method = CNR, Variance Estimator = Naive
    apply(t(sapply(simulation_list, function(xlist){
      sqrt(diag(xlist$cnr_out$naive_cov))
    })),2,mean) / apply(hatbeta_cnr_mat,2,sd),
    #Method = CNR, Variance Estimator = Naive BC
    apply(t(sapply(simulation_list, function(xlist){
      sqrt(diag(compute_naive_BC_cov_mat(cnr_out = xlist$cnr_out, pre_bs = xlist$pre_bs)))
    })),2,mean) / apply(hatbeta_cnr_mat,2,sd),
    #Method = CNR, Variance Estimator = Bootstrap
    apply(t(sapply(simulation_list,function(xlist){
      apply(xlist$second_bs$hatbeta_mat,2,sd)
    })),2,mean) / apply(hatbeta_cnr_mat,2,sd),
    #Method = CNR, Variance Estimator = Unadj-Bootstrap
    apply(t(sapply(simulation_list,function(xlist){
      apply(xlist$pre_bs$hatbeta_mat,2,sd)
    })),2,mean) / apply(hatbeta_cnr_mat,2,sd),
    #Method = CNR, Variance Estimator = NCC-Bootstrap
    apply(t(sapply(simulation_list,function(xlist){
      apply(xlist$NCC_bs$hatbeta_mat,2,sd)
    })),2,mean) / apply(hatbeta_cnr_mat,2,sd),
    #Method = BC-CNR, Variance Estimator = Naive
    apply(t(sapply(simulation_list, function(xlist){
      sqrt(diag(xlist$cnr_out$naive_cov))
    })),2,mean) / apply(hatbeta_bc_mat,2,sd),
    #Method = BC-CNR, Variance Estimator = Naive BC
    apply(t(sapply(simulation_list, function(xlist){
      sqrt(diag(compute_naive_BC_cov_mat(cnr_out = xlist$cnr_out, pre_bs = xlist$pre_bs)))
    })),2,mean) / apply(hatbeta_bc_mat,2,sd),
    #Method = BC-CNR, Variance Estimator = Bootstrap
    apply(t(sapply(simulation_list,function(xlist){
      apply(xlist$second_bs$hatbeta_mat,2,sd)
    })),2,mean) / apply(hatbeta_bc_mat,2,sd),
    #Method = BC-CNR, Variance Estimator = Unadj-Bootstrap
    apply(t(sapply(simulation_list,function(xlist){
      apply(xlist$pre_bs$hatbeta_mat,2,sd)
    })),2,mean) / apply(hatbeta_bc_mat,2,sd),
    #Method = BC-CNR, Variance Estimator = NCC-Bootstrap
    apply(t(sapply(simulation_list,function(xlist){
      apply(xlist$NCC_bs$hatbeta_mat,2,sd)
    })),2,mean) / apply(hatbeta_bc_mat,2,sd),
    #Method = 1-NMR, Variance Estimator = Naive
    apply(t(sapply(simulation_list, function(xlist){
      sqrt(diag(xlist$lnmr_L1$naive_cov))
    })),2,mean) / apply(hatbeta_1nmr_mat,2,sd),
    #Method = 3-NMR, Variance Estimator = Naive
    apply(t(sapply(simulation_list, function(xlist){
      sqrt(diag(xlist$lnmr_L3$naive_cov))
    })),2,mean) / apply(hatbeta_3nmr_mat,2,sd),
    #Method = 5-NMR, Variance Estimator = Naive
    apply(t(sapply(simulation_list, function(xlist){
      sqrt(diag(xlist$lnmr_L5$naive_cov))
    })),2,mean) / apply(hatbeta_5nmr_mat,2,sd)
  )
  
  #Extract point estimates of sigma2_rho, alpha_rho and tau_epsilon
  hattheta_err_cnr_mat = t(sapply(simulation_list,function(xlist){
    c(xlist$cnr_out$hatsigma2_rho,xlist$cnr_out$hatalpha_rho,xlist$cnr_out$hattau_epsilon)
  }))
  
  hattheta_err_bc_mat = t(sapply(simulation_list,function(xlist){
    c(exp(2 * log(xlist$cnr_out$hatsigma2_rho) - mean(log(xlist$pre_bs$hatsigma2_rho_vec), na.rm = T )   ),
      exp(2 * log(xlist$cnr_out$hatalpha_rho) - mean(log(xlist$pre_bs$hatalpha_rho_vec), na.rm = T )   ),
      exp(2 * log(xlist$cnr_out$hattau_epsilon) - mean(log(xlist$pre_bs$hattau_epsilon_vec), na.rm = T )   ))
  }))
  
  hattheta_err_1nmr_mat = t(sapply(simulation_list,function(xlist){
    c(xlist$lnmr_L1$hatsigma2_rho,xlist$lnmr_L1$hatalpha_rho,xlist$lnmr_L1$hattau_epsilon)
  }))
  
  hattheta_err_3nmr_mat = t(sapply(simulation_list,function(xlist){
    c(xlist$lnmr_L3$hatsigma2_rho,xlist$lnmr_L3$hatalpha_rho,xlist$lnmr_L3$hattau_epsilon)
  }))
  
  hattheta_err_5nmr_mat = t(sapply(simulation_list,function(xlist){
    c(xlist$lnmr_L5$hatsigma2_rho,xlist$lnmr_L5$hatalpha_rho,xlist$lnmr_L5$hattau_epsilon)
  }))
  
  bias_theta_err_mat = rbind(apply(hattheta_err_cnr_mat,2,mean) - c(true_sigma2_rho, true_alpha_rho, true_tau_epsilon),
                             apply(hattheta_err_bc_mat,2,mean) - c(true_sigma2_rho, true_alpha_rho, true_tau_epsilon),
                             apply(hattheta_err_1nmr_mat,2,mean) - c(true_sigma2_rho, true_alpha_rho, true_tau_epsilon),
                             apply(hattheta_err_3nmr_mat,2,mean) - c(true_sigma2_rho, true_alpha_rho, true_tau_epsilon),
                             apply(hattheta_err_5nmr_mat,2,mean) - c(true_sigma2_rho, true_alpha_rho, true_tau_epsilon))
  colnames(bias_theta_err_mat) = c('sigma2_rho','alpha_rho','tau_epsilon')
  
  var_theta_err_mat = rbind(apply(hattheta_err_cnr_mat,2,var),
                            apply(hattheta_err_bc_mat,2,var),
                            apply(hattheta_err_1nmr_mat,2,var),
                            apply(hattheta_err_3nmr_mat,2,var),
                            apply(hattheta_err_5nmr_mat,2,var))
  colnames(var_theta_err_mat) = c('sigma2_rho','alpha_rho','tau_epsilon')
  
  rmse_theta_err_mat = sqrt(bias_theta_err_mat^2 + var_theta_err_mat)
  colnames(rmse_theta_err_mat) = c('sigma2_rho','alpha_rho','tau_epsilon')
  
  #CI coverage
  
  #Method = CNR, CI Method = Naive
  hatbeta_cnr_naive_coverage = apply(t(sapply(simulation_list,function(xlist){
    (true_beta > xlist$cnr_out$hatbeta - qnorm(0.975) * sqrt(diag(xlist$cnr_out$naive_cov))) & 
      (true_beta < xlist$cnr_out$hatbeta + qnorm(0.975) * sqrt(diag(xlist$cnr_out$naive_cov)))
  })), 2, mean)
  
  #Method = CNR, CI Method = Naive BC
  hatbeta_cnr_naive_BC_coverage = apply(t(sapply(simulation_list,function(xlist){
    (true_beta > xlist$cnr_out$hatbeta - qnorm(0.975) * sqrt(diag(compute_naive_BC_cov_mat(cnr_out = xlist$cnr_out, pre_bs = xlist$pre_bs)))) & 
      (true_beta < xlist$cnr_out$hatbeta + qnorm(0.975) * sqrt(diag(compute_naive_BC_cov_mat(cnr_out = xlist$cnr_out, pre_bs = xlist$pre_bs))))
  })), 2, mean)
  
  #Method = CNR, CI Method = Bootstrap
  hatbeta_cnr_bs_quantile_coverage = apply(t(sapply(simulation_list,function(xlist){
    (true_beta > apply(xlist$second_bs$hatbeta_mat,2, function(x){quantile(x,probs = 0.025)}) ) & 
      (true_beta < apply(xlist$second_bs$hatbeta_mat,2, function(x){quantile(x,probs = 0.975)}) )
  })), 2, mean)
  
  #Method = CNR, CI Method = Unadj-Bootstrap
  hatbeta_cnr_unadj_bs_quantile_coverage = apply(t(sapply(simulation_list,function(xlist){
    (true_beta > apply(xlist$pre_bs$hatbeta_mat,2, function(x){quantile(x,probs = 0.025)}) ) & 
      (true_beta < apply(xlist$pre_bs$hatbeta_mat,2, function(x){quantile(x,probs = 0.975)}) )
  })), 2, mean)
  
  #Method = CNR, CI Method = NCC-Bootstrap
  hatbeta_cnr_NCC_bs_quantile_coverage = apply(t(sapply(simulation_list,function(xlist){
    (true_beta > apply(xlist$NCC_bs$hatbeta_mat,2, function(x){quantile(x,probs = 0.025)}) ) & 
      (true_beta < apply(xlist$NCC_bs$hatbeta_mat,2, function(x){quantile(x,probs = 0.975)}) )
  })), 2, mean)
  
  #Method = BC-CNR, CI Method = Naive
  hatbeta_bccnr_naive_coverage = apply(t(sapply(simulation_list,function(xlist){
    (true_beta > (2*xlist$cnr_out$hatbeta - xlist$pre_bs$bootstrap_hatbeta_mean) - qnorm(0.975) * sqrt(diag(xlist$cnr_out$naive_cov))) & 
      (true_beta < (2*xlist$cnr_out$hatbeta - xlist$pre_bs$bootstrap_hatbeta_mean) + qnorm(0.975) * sqrt(diag(xlist$cnr_out$naive_cov)))
  })), 2, mean)
  
  #Method = BC-CNR, CI Method = Naive BC
  hatbeta_bccnr_naive_BC_coverage = apply(t(sapply(simulation_list,function(xlist){
    (true_beta > (2*xlist$cnr_out$hatbeta - xlist$pre_bs$bootstrap_hatbeta_mean) - qnorm(0.975) * sqrt(diag(compute_naive_BC_cov_mat(cnr_out = xlist$cnr_out, pre_bs = xlist$pre_bs)))) & 
      (true_beta < (2*xlist$cnr_out$hatbeta - xlist$pre_bs$bootstrap_hatbeta_mean) + qnorm(0.975) * sqrt(diag(compute_naive_BC_cov_mat(cnr_out = xlist$cnr_out, pre_bs = xlist$pre_bs))))
  })), 2, mean)
  
  #Method = 1-NMR, CI Method = Naive
  hatbeta_1nmr_naive_coverage = apply(t(sapply(simulation_list,function(xlist){
    (true_beta > xlist$lnmr_L1$hatbeta - qnorm(0.975) * sqrt(diag(xlist$lnmr_L1$naive_cov))) & 
      (true_beta < xlist$lnmr_L1$hatbeta + qnorm(0.975) * sqrt(diag(xlist$lnmr_L1$naive_cov)))
  })), 2, mean)
  
  #Method = 3-NMR, CI Method = Naive
  hatbeta_3nmr_naive_coverage = apply(t(sapply(simulation_list,function(xlist){
    (true_beta > xlist$lnmr_L3$hatbeta - qnorm(0.975) * sqrt(diag(xlist$lnmr_L3$naive_cov))) & 
      (true_beta < xlist$lnmr_L3$hatbeta + qnorm(0.975) * sqrt(diag(xlist$lnmr_L3$naive_cov)))
  })), 2, mean)
  
  #Method = 5-NMR, CI Method = Naive
  hatbeta_5nmr_naive_coverage = apply(t(sapply(simulation_list,function(xlist){
    (true_beta > xlist$lnmr_L5$hatbeta - qnorm(0.975) * sqrt(diag(xlist$lnmr_L5$naive_cov))) & 
      (true_beta < xlist$lnmr_L5$hatbeta + qnorm(0.975) * sqrt(diag(xlist$lnmr_L5$naive_cov)))
  })), 2, mean)
  
  beta_CI_coverage_mat = rbind(hatbeta_cnr_naive_coverage, hatbeta_cnr_naive_BC_coverage, hatbeta_cnr_bs_quantile_coverage, hatbeta_cnr_unadj_bs_quantile_coverage, hatbeta_cnr_NCC_bs_quantile_coverage,
                               hatbeta_bccnr_naive_coverage, hatbeta_bccnr_naive_BC_coverage,
                               hatbeta_1nmr_naive_coverage, hatbeta_3nmr_naive_coverage, hatbeta_5nmr_naive_coverage)
  
  #Average CI Width
  #Method = CNR, CI Method = Naive
  hatbeta_cnr_naive_width = apply(t(sapply(simulation_list,function(xlist){
    2*qnorm(0.975) * sqrt(diag(xlist$cnr_out$naive_cov))
  })), 2, mean)
  
  #Method = CNR, CI Method = Naive BC
  hatbeta_cnr_naive_BC_width = apply(t(sapply(simulation_list,function(xlist){
    2*qnorm(0.975) * sqrt(diag(compute_naive_BC_cov_mat(cnr_out = xlist$cnr_out, pre_bs = xlist$pre_bs)))
  })), 2, mean)
  
  #Method = CNR, CI Method = Bootstrap
  hatbeta_cnr_bs_quantile_width = apply(t(sapply(simulation_list,function(xlist){
    apply(xlist$second_bs$hatbeta_mat,2, function(x){quantile(x,probs = 0.975)} ) -  apply(xlist$second_bs$hatbeta_mat,2, function(x){quantile(x,probs = 0.025)}) 
  })), 2, mean)
  
  #Method = CNR, CI Method = Unadj-Bootstrap
  hatbeta_cnr_unadj_bs_quantile_width = apply(t(sapply(simulation_list,function(xlist){
    apply(xlist$pre_bs$hatbeta_mat,2, function(x){quantile(x,probs = 0.975)} ) -  apply(xlist$pre_bs$hatbeta_mat,2, function(x){quantile(x,probs = 0.025)}) 
  })), 2, mean)
  
  #Method = CNR, CI Method = NCC-Bootstrap
  hatbeta_cnr_NCC_bs_quantile_width = apply(t(sapply(simulation_list,function(xlist){
    apply(xlist$NCC_bs$hatbeta_mat,2, function(x){quantile(x,probs = 0.975)} ) -  apply(xlist$NCC_bs$hatbeta_mat,2, function(x){quantile(x,probs = 0.025)}) 
  })), 2, mean)
  
  #Method = BC-CNR, CI Method = Naive
  hatbeta_bccnr_naive_width = apply(t(sapply(simulation_list,function(xlist){
    2*qnorm(0.975) * sqrt(diag(xlist$cnr_out$naive_cov))
  })), 2, mean)
  
  #Method = BC-CNR, CI Method = Naive BC
  hatbeta_bccnr_naive_BC_width = apply(t(sapply(simulation_list,function(xlist){
    2*qnorm(0.975) * sqrt(diag(compute_naive_BC_cov_mat(cnr_out = xlist$cnr_out, pre_bs = xlist$pre_bs)))
  })), 2, mean)
  
  #Method = 1-NMR, CI Method = Naive
  hatbeta_1nmr_naive_width = apply(t(sapply(simulation_list,function(xlist){
    2*qnorm(0.975) * sqrt(diag(xlist$lnmr_L1$naive_cov))
  })), 2, mean)
  
  
  #Method = 3-NMR, CI Method = Naive
  hatbeta_3nmr_naive_width = apply(t(sapply(simulation_list,function(xlist){
    2*qnorm(0.975) * sqrt(diag(xlist$lnmr_L3$naive_cov))
  })), 2, mean)
  
  #Method = 5-NMR, CI Method = Naive
  hatbeta_5nmr_naive_width = apply(t(sapply(simulation_list,function(xlist){
    2*qnorm(0.975) * sqrt(diag(xlist$lnmr_L5$naive_cov))
  })), 2, mean)
  
  beta_CI_average_width_mat = rbind(hatbeta_cnr_naive_width, hatbeta_cnr_naive_BC_width, hatbeta_cnr_bs_quantile_width, hatbeta_cnr_unadj_bs_quantile_width, hatbeta_cnr_NCC_bs_quantile_width,
                                    hatbeta_bccnr_naive_width, hatbeta_bccnr_naive_BC_width,
                                    hatbeta_1nmr_naive_width, hatbeta_3nmr_naive_width, hatbeta_5nmr_naive_width)
  
  bias_beta_df = data.frame(Method = c('CNR','BC-CNR','1-NMR','3-NMR','5-NMR'),
                            bias_beta_mat[,-1])
  
  rmse_beta_df = data.frame(Method = c('CNR','BC-CNR','1-NMR','3-NMR','5-NMR'),
                            rmse_beta_mat[,-1])
  
  ase_esd_ratio_df = data.frame(Method = c(rep('CNR',5), rep('BC-CNR',5), '1-NMR', '3-NMR', '5-NMR'),
                                Variance_estimator = c( rep(c('Naive', 'Naive BC', 'Bootstrap', 'Unadj-Bootstrap', 'NCC-Bootstrap'),2), rep('Naive',3) ),
                                ase_esd_ratio_mat[,-1])
  
  bias_theta_err_df = data.frame(Method = c('CNR','BC-CNR','1-NMR','3-NMR','5-NMR'),
                                 bias_theta_err_mat)
  
  rmse_theta_err_df = data.frame(Method = c('CNR','BC-CNR','1-NMR','3-NMR','5-NMR'),
                                 rmse_theta_err_mat)
  
  beta_CI_coverage_df = data.frame(Method = c(rep('CNR',5), rep('BC-CNR',2), '1-NMR', '3-NMR', '5-NMR'),
                                   CI_method = c('Naive', 'Naive BC', 'Bootstrap', 'Unadj-Bootstrap', 'NCC-Bootstrap', 'Naive', 'Naive BC', rep('Naive',3)),
                                   beta_CI_coverage_mat[,-1])
  
  beta_CI_average_width_df = data.frame(Method = c(rep('CNR',5), rep('BC-CNR',2), '1-NMR', '3-NMR', '5-NMR'),
                                        CI_method = c('Naive', 'Naive BC', 'Bootstrap', 'Unadj-Bootstrap', 'NCC-Bootstrap', 'Naive', 'Naive BC', rep('Naive',3)),
                                        beta_CI_average_width_mat[,-1])
  
  return(list(
    bias_beta_df = bias_beta_df, rmse_beta_df = rmse_beta_df, ase_esd_ratio_df = ase_esd_ratio_df, 
    bias_theta_err_df = bias_theta_err_df, rmse_theta_err_df = rmse_theta_err_df,
    beta_CI_coverage_df = beta_CI_coverage_df, beta_CI_average_width_df = beta_CI_average_width_df
  ))
  
}

#' @title 
#' Compute the naive BC covariance matrix of CNR mean regression coefficient estimates.
#' 
#' @description 
#' Compute the naive BC covariance matrix of CNR mean regression coefficient estimates based on bias-corrected version of estimated covariance matrix of spatial random effects plus residuals.
#' 
#' @param cnr_out A list returned by the \code{cnr} function.
#' @param pre_bs A list returned by the \code{bootstrap_cnr} function, representing the preliminary bootstrap results.
#' 
#' @return the naive BC covariance matrix of CNR mean regression coefficient estimates based on bias-corrected version of estimated covariance matrix of spatial random effects plus residuals

compute_naive_BC_cov_mat = function(cnr_out, pre_bs){
  m = nrow(cnr_out$data_x); n = nrow(cnr_out$data_y)
  
  #Bias-correction
  bc_hatsigma2_rho = exp(2 * log(cnr_out$hatsigma2_rho) - mean(log(pre_bs$hatsigma2_rho_vec), na.rm = T )   )
  bc_hatnu_rho = exp(2 * log(cnr_out$hatnu_rho) - mean(log(pre_bs$hatnu_rho_vec), na.rm = T )   )
  bc_hatalpha_rho = exp(2 * log(cnr_out$hatalpha_rho) - mean(log(pre_bs$hatalpha_rho_vec), na.rm = T )   )
  bc_hattau_epsilon = exp(2 * log(cnr_out$hattau_epsilon) - mean(log(pre_bs$hattau_epsilon_vec), na.rm = T )   )
  
  bc_hatSigma_err = as.numeric(bc_hatsigma2_rho) * vec2symMat(geoR::matern(u = cnr_out$dist_vec , kappa = bc_hatnu_rho, phi = (1/bc_hatalpha_rho)), diag = F ) + bc_hattau_epsilon * diag(m+n)
  bc_hatSigma_err = bc_hatSigma_err[(m+1):(m+n),(m+1):(m+n)]
  bc_naive_cov = solve( t(cnr_out$step3_fitme$X.pv) %*% solve(bc_hatSigma_err) %*% cnr_out$step3_fitme$X.pv )
  return(bc_naive_cov)
}
