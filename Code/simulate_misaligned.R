#' @title 
#' Simulating spatially misaligned data
#' 
#' @description 
#' Generate spatially misaligned response and covariates, where the covariates follow a multivariate Gaussian distribution with a generalized Kronecker product covariance structure consisting of marginal Mat\'{e}rn covariance matrix and an unstructured cross-correlation matrix. 
#' The response is generated from a spatial linear mixed model, where the mean structure is some function of the covariates and the spatial random effects are also assumed to have a Mat\'{e}rn covariance matrix. 
#' 
#' @details 
#' This function can be used to generate spatially misaligned data for either simulation (\code{bootstrap = FALSE}) or bootstrap purpose (\code{bootstrap = TRUE}). 
#' For simulation purpose, \code{spline_basis_list} should be specified when the mean structure of the response involves splines of the covariates, otherwise \code{y_mean_model_formula} should be specified instead,
#' noting that this formula *should not* include the response on the LHS and the spatial random effect on the RHS. For example, if \code{x_name = c('X1','X2')} and the mean structure is a linear combination of 'X1' and 'X2', then \code{y_mean_model_formula = as.formula('~ X1 + X2')}.
#' For bootstrap purpose, the function requires pre-computed lower Cholesky factors \code{tchol_Sigma_mat} and \code{tchol_Sigma_err}, together with the \code{spaMM::fitme} model object returned by \code{cnr_out$step3_fitme} where \code{cnr_out} is the returned list from \code{cnr} function.
#' 
#' @param nu_vec_x A vector of Mat\'{e}rn smoothness parameters for the covariates. This has to be specified when \code{bootstrap_use = FALSE}.
#' @param alpha_vec_x A vector of Mat\'{e}rn range parameters for the covariates. This has to be specified when \code{bootstrap_use = FALSE}.
#' @param sigma2_vec_x A vector of Mat\'{e}rn variance parameters for the covariates. This has to be specified when \code{bootstrap_use = FALSE}.
#' @param tau_vec_x A vector of nugget parameters for the covariates. This has to be specified when \code{bootstrap_use = FALSE}.
#' @param mu_vec_x A vector of mean parameters for the covariates.
#' @param R_mat A cross-correlation matrix for the covariates. This has to be specified when \code{bootstrap_use = FALSE}.
#' @param beta_vec A vector of mean regression coefficients.
#' @param sigma2_rho Varianc parameter of spatial random effect. This has to be specified when \code{bootstrap_use = FALSE}.
#' @param nu_rho Mat\'{e}rn smoothness parameter for the spatial random effect. This has to be specified when \code{bootstrap_use = FALSE}.
#' @param alpha_rho Mat\'{e}rn range parameter for the spatial random effect. This has to be specified when \code{bootstrap_use = FALSE}.
#' @param tau_epsilon Residual variance parameter, or equivalently, the nugget parameter for the response. This has to be specified when \code{bootstrap_use = FALSE}.
#' @param y_loc A data frame describing the response locations. It should consist of two columns named "lon" and "lat" denoting longitude and latitude, respectively, when \code{gc_dist = TRUE}.
#' @param x_loc A data frame describing the misaligned covariate locations. It should consist of two columns named "lon" and "lat" denoting longitude and latitude, respectively, when \code{gc_dist = TRUE}.
#' @param y_name A character string representing the name of the response. 
#' @param x_name A character vector containing the names of the covariates.
#' @param y_mean_model_formula A model formula indicating the mean structure of the response as function of the covariates, the details of this is given under 'Details'. This or \code{spline_basis_list} has to be specified when \code{bootstrap_use = FALSE}.
#' @param spline_basis_list A list with \code{length(x_name)} number of elements, where each element consists of the spline basis matrix. This is used when splines (of the covariates) are used to model the mean structure of the response. This or \code{y_mean_model_formula} has to be specified when \code{bootstrap_use = FALSE}.
#' @param cnr_step3_model_fit A \code{spaMM::fitme} model object returned by \code{cnr_out$step3_fitme} where \code{cnr_out} is the returned list from \code{cnr} function. This has to be specified when \code{bootstrap_use = TRUE}.
#' @param gc_dist Logical. If \code{TRUE} great-circle distance (in kilometers) is used in the Mat\'{e}rn covariance matrices, otherwise Euclidean distance is used. This has to be specified when \code{bootstrap_use = FALSE}.
#' @param bootstrap_use Logical. If \code{TRUE} it indicates that the function is being used in Step (i) of the bootstrap procedure, and users should supply pre-computed lower Cholesky factors \code{tchol_Sigma_mat} and \code{tchol_Sigma_err} to save computational time. 
#' @param tchol_Sigma_mat A pre-computed lower Cholesky factor of the estimated joint covariance matrix of all covariates, which can be obtained from \code{compute_tchol_hatSigma_mat}. This has to be specified when \code{bootstrap_use = TRUE}.
#' @param tchol_Sigma_err A pre-computed lower Cholesky factor of the estimated covariance matrix of the spatial random effect plus the residuals, which can be obtained from \code{compute_tchol_hatSigma_err}. This has to be specified when \code{bootstrap_use = TRUE}.
#'
#' @return A list with the elements being input arguments of the function call and the following additional elements:
#' \item{data_x:}{A data frame consisting of the covariates at the misaligned covariate locations.}
#' \item{data_y:}{A data frame consisting of the response at the response locations.}
#' \item{loc_name:}{A character vector containing the names of the coordinates. When \code{gc_dist = TRUE}, it should be equal to \code{c('lon','lat')}.}



simulate_misaligned = function(nu_vec_x = NULL, alpha_vec_x = NULL, sigma2_vec_x = NULL, tau_vec_x = NULL, mu_vec_x, R_mat = NULL, beta_vec, sigma2_rho = NULL, nu_rho = NULL, alpha_rho = NULL, tau_epsilon = NULL, y_loc, x_loc, y_name, x_name, y_mean_model_formula = NULL , spline_basis_list = NULL, cnr_step3_model_fit = NULL, gc_dist = NULL,
                               bootstrap_use = FALSE, tchol_Sigma_mat = NULL, tchol_Sigma_err = NULL){
  m = dim(x_loc)[1]; n = dim(y_loc)[1]; K = length(x_name)
  if(!bootstrap_use){
    all_loc = rbind(x_loc,y_loc)
    if(gc_dist){
      all_loc_sf = st_as_sf(all_loc, coords = c("lon", "lat"), crs = "WGS84", agr = "constant")
      dist_mat =  st_distance(all_loc_sf)
      dist_mat = dist_mat / 1000
    }else{
      dist_mat = as.matrix(dist(all_loc))
      dist_mat_x = as.matrix(dist(x_loc))
    }
    dist_vec = dist_mat[lower.tri(dist_mat)]
    dist_vec = as.numeric(dist_vec)
    
    Sigma_err = sigma2_rho * vec2symMat(geoR::matern(u = dist_vec , kappa = nu_rho, phi = (1/alpha_rho)), diag = F ) + tau_epsilon * diag(m+n)
    tchol_Sigma_err = t(chol(Sigma_err[(m+1):(m+n),(m+1):(m+n)]))
    
    L_k_list = lapply(1:K, function(k){
      t(chol(sigma2_vec_x[k] * vec2symMat(geoR::matern(u = dist_vec, kappa = nu_vec_x[k], phi = (1/(alpha_vec_x[k] ))), diag = F ) + tau_vec_x[k] * diag(m+n)))
    })
    Sigma_mat = as.matrix(bdiag(L_k_list)  %*% (R_mat %x% diag(m+n) ) %*%  t( bdiag(L_k_list) ) )
    tchol_Sigma_mat = t(chol(Sigma_mat))
  }
  
  #Generating covariates where x ~ N (mu_vec_x \otimes 1_{m+n}, Sigma_mat)
  x = matrix( (mu_vec_x %x% rep(1,m+n)), ncol = 1) +  ( tchol_Sigma_mat %*% rnorm((m+n) *K ) )
  x_tildeS = data.frame((matrix(x, nrow = m+n, ncol = K, byrow = FALSE, dimnames = list(c(),x_name)))[1:m,])
  x_S = data.frame((matrix(x, nrow = m+n, ncol = K, byrow = FALSE, dimnames = list(c(),x_name)))[(m+1):(m+n),])
  
  #Need y_mean_model_formula to have some form like ' ~ BLABLABLA ' without the random effect
  # or need the basis matrices for each of the covariate - which is used when this function is used in step (i) of the bootstrap procedure
  #spline_basis_list should have the same name as the columns of x_S
  if(bootstrap_use){
    y = as.vector(predict(object = cnr_step3_model_fit, newdata = x_S, re.form = NA)) + tchol_Sigma_err %*% rnorm(n)
  }else{
    if(is.null(y_mean_model_formula)){
      temp_X = do.call("cbind",lapply(1:K, function(k){
        predict(object = spline_basis_list[[ x_name[k] ]], newx = x_S[,x_name[k] ])
      }))
      y = cbind(1,temp_X) %*% beta_vec + tchol_Sigma_err %*% rnorm(n)
    }else{
      y = model.matrix(y_mean_model_formula, data = x_S) %*% beta_vec + tchol_Sigma_err %*% rnorm(n)
    }
  }
  
  data_x = data.frame(x_tildeS, x_loc)
  data_y = data.frame(y, y_loc)
  colnames(data_y) = c(y_name, colnames(y_loc))
  
  return(list(
    #Input arguments
    nu_vec_x = nu_vec_x, alpha_vec_x = alpha_vec_x, sigma2_vec_x = sigma2_vec_x, tau_vec_x = tau_vec_x, mu_vec_x = mu_vec_x, R_mat = R_mat, beta_vec = beta_vec, sigma2_rho = sigma2_rho, nu_rho = nu_rho, alpha_rho = alpha_rho, tau_epsilon = tau_epsilon, y_loc = y_loc, x_loc = x_loc, y_mean_model_formula = y_mean_model_formula, y_name = y_name, x_name = x_name, spline_basis_list = spline_basis_list ,gc_dist = gc_dist,
    #Additional stuffs
    data_x = data_x, data_y = data_y, loc_name = colnames(y_loc)
  ))
}