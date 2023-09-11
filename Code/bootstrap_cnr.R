#' @title 
#' Parametric bootstrap for cokrig-and-regress
#' 
#' @description 
#' Perform parametric bootstrap for bias-correction or uncertainty quantification for cokrig-and-regress
#' 
#' @details 
#' This function perform 
#' T_bootstrap, bootstrap_CI_level, cnr_out, bootstrap_type, pre_bs = NULL
#' @param T_bootstrap Number of bootstrap datasets.
#' @param bootstrap_CI_level Confidence level for constructing the bootstrap percentile intervals of the mean regression coefficients and conditional smoothers (if \code{!is.null(cnr_out$smooth_list)}) with the default being 0.95.
#' @param cnr_out A list returned by the \code{cnr} function.
#' @param bootstrap_type Can be one of 'Prelim', 'Second', and 'NCC', representing the preliminary, second and non-cross-correlated bootstraps, respetively.
#' @param pre_bs A list returned by this function, representing the preliminary bootstrap results. This only has to be specified when \code{bootstrap_type %in% c('Second','NCC')}.
#' 
#' @return A list with the following elements:
#' \item{hatnu_mat_x:}{A matrix with each row representing the Mat\'{e}rn smoothness parameters of the covariates that are assumed known for each bootstrap sample.}
#' \item{hatalpha_mat_x:}{A matrix with each row representing the estimated Mat\'{e}rn range parameters of the covariates for each bootstrap sample.}
#' \item{hatsigma2_mat_x:}{A matrix with each row representing the estimated Mat\'{e}rn variance parameter of the covariates for each bootstrap sample.}
#' \item{hattau_mat_x:}{A matrix with each row representing the estimated nugget parameters of the covariates for each bootstrap sample.}
#' \item{hatmu_mat_x:}{A matrix with each row representing the estimated mean parameters of the covariates for each bootstrap sample.}
#' \item{hatR_mat_x:}{A matrix with each row representing the lower triangular part of the estimated cross-correlation matrix of the covariates for each bootstrap sample.}
#' \item{hatbeta_mat_x:}{A matrix with each row representing the estimated mean regression coefficients for each bootstrap sample.}
#' \item{hatsigma2_rho_vec:}{A vector with each element representing the estimated Mat\'{e}rn variance parameter of the spatial random effect for each bootstrap sample.}
#' \item{hatnu_rho_vec:}{A vector with each element representing the Mat\'{e}rn smoothness parameter of the spatial random effect that is assumed known for each bootstrap sample.}
#' \item{hatalpha_rho_vec:}{A vector with each element representing the estimated Mat\'{e}rn range parameter of the spatial random effect for each bootstrap sample.}
#' \item{hattau_epsilon_vec:}{A vector with each element representing the estimated residual variance parameter for each bootstrap sample.}
#' \item{step1_warning_list:}{A list with each element representing another list of warnings for events that may have occurred during the \code{spaMM::fitme} fit in Step (ii) of each bootstrap sample.}
#' \item{step1_fitme_time_total_vec:}{A vector with each element representing the total time taken for the \code{spaMM::fitme} fit in Step (ii) of each bootstrap sample.}
#' \item{step3_warning_list:}{A list with each element representing another list of warnings for events that may have occurred during the \code{spaMM::fitme} fit in Step (iv) of each bootstrap sample.}
#' \item{step3_fitme_time_vec:}{A vector with each element representing the total time taken for the \code{spaMM::fitme} fit in Step (iv) of each bootstrap sample.}
#' \item{bootstrap_hatbeta_mean:}{A vector of bootstrap mean for the mean regression coefficients.}
#' \item{bootstrap_hatsigma2_rho_mean:}{Bootstrap mean for the Mat\'{e}rn variance parameter of the spatial random effect.}
#' \item{bootstrap_hatnu_rho_mean:}{Bootstrap mean for the Mat\'{e}rn smoothness parameter of the spatial random effect that is assumed known.}
#' \item{bootstrap_hatalpha_rho_mean:}{Bootstrap mean for the Mat\'{e}rn range parameter of the spatial random effect.}
#' \item{bootstrap_hattau_epsilon_mean:}{Bootstrap mean for the residual variance parameter.}
#' \item{bootstrap_hatbeta_se:}{A vector of bootstrap standard deviation for the mean regression coefficients.}
#' \item{bootstrap_percentile_CI_beta:}{A matrix representing the bootstrap percentile confidence intervals for the mean regression coefficients, where the first row is the lower confidence limit and the second row is the upper confidence limit.}
#' \item{bootstrap_smoother:}{If \code{!is.null(cnr_out$smooth_list)}, then it is a list with \code{length(cnr_out$x_name)} number of elements, and each element of the list is a data frame consisting of \code{3 + length(cnr_out$x_name)} columns. In particular, these data frames are similar to the elements of \code{cnr_out$smoother}, with the difference being the first three columns in the elements of \code{bootstrap_smoother}. As opposed to \code{cnr_out$smoother}, the first column represents the bootstrap means of the conditional smoothers, while the second and third column represent the lower and upper confidence limits of the bootstrap percentile interval of the conditional smoothers. If \code{is.null(cnr_out$smooth_list)}, then it is \code{NULL}.}
#' 
#' 

bootstrap_cnr = function(T_bootstrap, bootstrap_CI_level = 0.95, cnr_out, bootstrap_type, pre_bs = NULL){
  m = dim(cnr_out$x_loc)[1]; n = dim(cnr_out$y_loc)[1]; K = length(cnr_out$x_name)
  
  
  #Initializing matrices/vectors/lists/arrays to store bootstrap results
  hatnu_mat_x = matrix(NA, nrow = T_bootstrap, ncol = K)
  hatalpha_mat_x = matrix(NA, nrow = T_bootstrap, ncol = K)
  hatsigma2_mat_x = matrix(NA, nrow = T_bootstrap, ncol = K)
  hattau_mat_x = matrix(NA, nrow = T_bootstrap, ncol = K)
  hatmu_mat_x = matrix(NA, nrow = T_bootstrap, ncol = K)
  hatR_mat_x = matrix(NA, nrow = T_bootstrap, ncol = (K*(K-1))/2) 
  
  step1_warning_list = vector("list", T_bootstrap)
  step1_fitme_time_total_vec = rep(NA, T_bootstrap)
  step3_warning_list = vector("list", T_bootstrap)
  step3_fitme_time_vec = rep(NA, T_bootstrap)
  
  
  hatbeta_mat = matrix(NA, nrow = T_bootstrap, ncol = length(cnr_out$hatbeta))
  hatsigma2_rho_vec = rep(NA, T_bootstrap)
  hatnu_rho_vec = rep(NA,T_bootstrap)
  hatalpha_rho_vec = rep(NA,T_bootstrap)
  hattau_epsilon_vec = rep(NA,T_bootstrap)
  
  if(!is.null(cnr_out$smooth_list)){
    smooth_matrix_array = array(NA, dim = c(nrow(cnr_out$smoother[[1]]),K,T_bootstrap), dimnames = list(c(),cnr_out$x_name, c() ) )
  }else{
    smooth_matrix_array = NULL
  }
  
  if(bootstrap_type %in% c('Prelim', 'Second')){
    tchol_hatSigma_mat = compute_tchol_hatSigma_mat(cnr_out = cnr_out, NCC = FALSE)
  }else{
    tchol_hatSigma_mat = compute_tchol_hatSigma_mat(cnr_out = cnr_out, NCC = TRUE)
  }
  
  if(bootstrap_type %in% c('Second', 'NCC')){
    tchol_hatSigma_err = compute_tchol_hatSigma_err(cnr_out = cnr_out, bc = TRUE, pre_bs = pre_bs)
  }else{
    tchol_hatSigma_err = compute_tchol_hatSigma_err(cnr_out = cnr_out, bc = FALSE, pre_bs = NULL)
  }
  
  for(t in 1:T_bootstrap){
    tryCatch({
      #Bootstrap Step (i)  - Generating bootstrap datasets
      bootstrap_data_t = simulate_misaligned(mu_vec_x = cnr_out$hatmu_vec_x, beta_vec = cnr_out$hatbeta, y_loc = cnr_out$data_y[, cnr_out$loc_name], x_loc = cnr_out$data_x[, cnr_out$loc_name], y_name = cnr_out$y_name, x_name = cnr_out$x_name, y_mean_model_formula = NULL, spline_basis_list = NULL, cnr_step3_model_fit = cnr_out$step3_fitme,
                                             bootstrap_use = TRUE, 
                                             tchol_Sigma_mat = tchol_hatSigma_mat, 
                                             tchol_Sigma_err = tchol_hatSigma_err)
      
      #Bootstrap Step (ii) - (iv) - Estimating hattheta_x, computing cokriging predictions, and finally fitting spatial linear mixed model
      bootstrap_result_t = cnr(data_x = bootstrap_data_t$data_x, 
                               data_y = bootstrap_data_t$data_y,  
                               gc_dist = cnr_out$gc_dist, 
                               nu_vec_x = cnr_out$nu_vec_x, 
                               dist_vec = cnr_out$dist_vec, 
                               y_name = cnr_out$y_name, 
                               x_name = cnr_out$x_name, 
                               loc_name = bootstrap_data_t$loc_name,
                               slmm_formula = cnr_out$slmm_formula, 
                               nu_rho = cnr_out$nu_rho, 
                               smooth_list = cnr_out$smooth_list_out, discretized_num = NULL, smooth_CI_level = NULL)
      
      hatnu_mat_x[t,] = bootstrap_result_t$hatnu_vec_x
      hatalpha_mat_x[t,] = bootstrap_result_t$hatalpha_vec_x
      hatsigma2_mat_x[t,] = bootstrap_result_t$hatsigma2_vec_x
      hattau_mat_x[t,] = bootstrap_result_t$hattau_vec_x
      hatmu_mat_x[t,] = bootstrap_result_t$hatmu_vec_x
      hatR_mat_x = bootstrap_result_t$hatR[lower.tri(bootstrap_result_t$hatR)]
      
      step1_warning_list[[t]] = bootstrap_result_t$step1_warnings
      step1_fitme_time_total_vec[t] = bootstrap_result_t$step1_fitme_time_total
      step3_warning_list[[t]] = bootstrap_result_t$step3_warnings
      step3_fitme_time_vec[t]= bootstrap_result_t$step3_fitme_time
      
      
      hatbeta_mat[t,] = bootstrap_result_t$hatbeta
      hatsigma2_rho_vec[t] = bootstrap_result_t$hatsigma2_rho
      hatnu_rho_vec[t] = bootstrap_result_t$hatnu_rho
      hatalpha_rho_vec[t] = bootstrap_result_t$hatalpha_rho
      hattau_epsilon_vec[t] = bootstrap_result_t$hattau_epsilon
      
      
      
      #Compute conditional smoothers if required
      if(!is.null(cnr_out$smooth_list)){
        #Only store the point estimates for the smoothers
        smooth_matrix_array[,,t] = sapply(bootstrap_result_t$smoother, function(xlist){
          xlist[,1]
        })
      }
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  #Compute bootstrap mean of bootstrap samples 
  bootstrap_hatbeta_mean = apply(hatbeta_mat,2,function(x){mean(x, na.rm = T)})
  bootstrap_hatsigma2_rho_mean = mean(hatsigma2_rho_vec, na.rm = T)
  bootstrap_hatnu_rho_mean = mean(hatnu_rho_vec, na.rm = T)
  bootstrap_hatalpha_rho_mean = mean(hatalpha_rho_vec, na.rm = T)
  bootstrap_hattau_epsilon_mean = mean(hattau_epsilon_vec, na.rm = T)
  
  #Compute bootstrap SE of bootstrap samples for hatbeta
  bootstrap_hatbeta_se = apply(hatbeta_mat,2,function(x){sd(x, na.rm = T)})
  
  #Compute bootstrap percentile intervals for beta based on the bootstrap_CI_level
  bootstrap_percentile_CI_beta = apply(hatbeta_mat,2,function(x){quantile(x, probs = c( (1-bootstrap_CI_level)/2  , 1 - (1-bootstrap_CI_level)/2), na.rm = T)})
  
  #Compute bootstrap mean and bootstrap percentile intervals for smoothers based on the bootstrap_CI_level
  if(!is.null(cnr_out$smooth_list)){
    bootstrap_smoother = lapply(1:K, function(k){
      df_out = data.frame( apply(smooth_matrix_array[,k,], 1,function(x){mean(x,na.rm = T)}),
                           apply(smooth_matrix_array[,k,], 1,function(x){quantile(x, probs = (1-bootstrap_CI_level)/2,na.rm = T)}),
                           apply(smooth_matrix_array[,k,], 1,function(x){quantile(x, probs = 1 - (1-bootstrap_CI_level)/2,na.rm = T)}),
                           cnr_out$smoother[[k]][, 4:ncol(cnr_out$smoother[[k]])] )
      colnames(df_out) = c(cnr_out$y_name, paste0(cnr_out$y_name,'_lower'), paste0(cnr_out$y_name,'_upper'), cnr_out$x_name)
      df_out
    })
  }else{
    bootstrap_smoother = NULL
  }
  
  return(list(
    hatnu_mat_x = hatnu_mat_x, hatalpha_mat_x = hatalpha_mat_x, hatsigma2_mat_x = hatsigma2_mat_x, hattau_mat_x = hattau_mat_x, hatmu_mat_x = hatmu_mat_x, hatR_mat_x = hatR_mat_x,
    hatbeta_mat = hatbeta_mat, hatsigma2_rho_vec = hatsigma2_rho_vec, hatnu_rho_vec = hatnu_rho_vec, hatalpha_rho_vec = hatalpha_rho_vec, hattau_epsilon_vec= hattau_epsilon_vec,
    step1_warning_list = step1_warning_list, step1_fitme_time_total_vec = step1_fitme_time_total_vec, step3_warning_list = step3_warning_list, step3_fitme_time_vec = step3_fitme_time_vec,
    bootstrap_hatbeta_mean = bootstrap_hatbeta_mean, bootstrap_hatsigma2_rho_mean = bootstrap_hatsigma2_rho_mean, bootstrap_hatnu_rho_mean = bootstrap_hatnu_rho_mean, bootstrap_hatalpha_rho_mean = bootstrap_hatalpha_rho_mean, bootstrap_hattau_epsilon_mean = bootstrap_hattau_epsilon_mean, 
    bootstrap_hatbeta_se = bootstrap_hatbeta_se, bootstrap_percentile_CI_beta = bootstrap_percentile_CI_beta, bootstrap_smoother = bootstrap_smoother
  ))
  
}

#' @title 
#' Pre-compute lower Cholesky factor of the CNR estimated covariance matrix for all covariates
#' 
#' @description 
#' Pre-compute lower Cholesky factor of the CNR estimated covariance matrix for all covariates, which can then be used to save computational effort for bootstrap.
#' 
#' @param cnr_out A list returned by the \code{cnr} function.
#' @param NCC Logical. If \code{NCC = TRUE}, then the estimated cross-correlation matrix is set to be an identity matrix, i.e., ignoring the estimated cross-correlation between the covariates, before computing the lower Cholesky factor. This is set to be \code{TRUE} only when pre-computing the lower Cholesky factor for non-cross-correlated bootstrap.
#' 
#' @return the lower Cholesky factor of the CNR estimated covariance matrix of all covariates.

compute_tchol_hatSigma_mat = function(cnr_out, NCC = FALSE){
  K = length(cnr_out$hatnu_vec_x)
  m_n = dim(vec2symMat(cnr_out$dist_vec, diag = F))[1]
  hatL_k_list = lapply(1:K, function(k){
    t(chol(cnr_out$hatsigma2_vec_x[k] * vec2symMat(geoR::matern(u = cnr_out$dist_vec, kappa = cnr_out$hatnu_vec_x[k], phi = (1/(cnr_out$hatalpha_vec_x[k] ))), diag = F ) + cnr_out$hattau_vec_x[k] * diag(m_n)))
  })
  if(!NCC){
    # hatR = cnr_out$hatR
    # diag(hatR) = 1
    # hatSigma_mat = as.matrix(bdiag(hatL_k_list)  %*% (hatR %x% diag(m_n) ) %*%  t( bdiag(hatL_k_list) ) )
    hatSigma_mat = as.matrix(bdiag(hatL_k_list)  %*% (cnr_out$hatR %x% diag(m_n) ) %*%  t( bdiag(hatL_k_list) ) )
  }else{
    hatSigma_mat = as.matrix(bdiag(hatL_k_list)  %*% (diag(K) %x% diag(m_n) ) %*%  t( bdiag(hatL_k_list) ) )
  }
  return(t(chol(hatSigma_mat)))
}

#' @title 
#' Pre-compute lower Cholesky factor of the CNR estimated covariance matrix for spatial random effects plus residuals.
#' 
#' @description 
#' Pre-compute lower Cholesky factor of the CNR estimated covariance matrix for spatial random effects plus residuals., which can then be used to save computational effort for bootstrap.
#' 
#' @param cnr_out A list returned by the \code{cnr} function.
#' @param bc Logical. If \code{bc = TRUE}, then the CNR estimates of the spatial covariance parameters are bias-corrected at the log scale using the preliminary bootstrap samples from \code{pre_bs}, before being used to compute the lower Cholesky factor. This is set to be \code{TRUE} only when pre-computing the lower Cholesky factor for second or non-cross-correlated bootstrap.
#' @param pre_bs A list returned by the \code{bootstrap_cnr} function, representing the preliminary bootstrap results. This only has to be specified when \code{bc = TRUE}.
#' 
#' @return the lower Cholesky factor of the CNR estimated covariance matrix for spatial random effects plus residuals.


compute_tchol_hatSigma_err = function(cnr_out, bc = FALSE, pre_bs = NULL){
  m = nrow(cnr_out$data_x); n = nrow(cnr_out$data_y)
  if(!bc){
    tchol_hatSigma_err = t(chol(cnr_out$hatSigma_err))
  }else{
    bc_hatsigma2_rho = exp(2 * log(cnr_out$hatsigma2_rho) - mean(log(pre_bs$hatsigma2_rho_vec), na.rm = T )   )
    bc_hatnu_rho = exp(2 * log(cnr_out$hatnu_rho) - mean(log(pre_bs$hatnu_rho_vec), na.rm = T )   )
    bc_hatalpha_rho = exp(2 * log(cnr_out$hatalpha_rho) - mean(log(pre_bs$hatalpha_rho_vec), na.rm = T )   )
    bc_hattau_epsilon = exp(2 * log(cnr_out$hattau_epsilon) - mean(log(pre_bs$hattau_epsilon_vec), na.rm = T )   )
    
    bc_hatSigma_err = bc_hatsigma2_rho * vec2symMat(geoR::matern(u = cnr_out$dist_vec, kappa = bc_hatnu_rho, phi = (1/bc_hatalpha_rho) ), diag = F )  + bc_hattau_epsilon * diag(m+n)
    bc_hatSigma_err = bc_hatSigma_err[(m+1):(m+n),(m+1):(m+n)]
    tchol_hatSigma_err = t(chol(bc_hatSigma_err))
  }
  return(tchol_hatSigma_err)
}