#' @title 
#' L-nearest-matching-and-regress on spatially misaligned data.
#' 
#' @description 
#' Perform l-nearest-matching-and-regress on spatially misaligned data.
#' 
#' @details 
#' This function performs l-nearest-matching-and-regress (l-NMR) on spatially misaligned data, where the continuous response is spatially misaligned with multiple covariates. 
#' For \code{dist_vec}, the relevant distance matrix \code{dist_mat} should be constructed by first combining \code{x_loc} and \code{y_loc} (with \code{x_loc} ordered first) and then computing the distance matrix based on all locations, where again the distance measure should be great-circle distance if \code{gc_dist = TRUE} or Euclidean distance if \code{gc_dist = FALSE}.
#' The model formula \code{slmm_formula} should be a linear model formula that can be handled by \code{spaMM::fitme} fitting function with a spatial random effect term that is assumed to have a Mat\'{e}rn covariance structure. For instance, \code{slmm_formula = as.formula('y ~ X1 + X2 + Matern(1| lon + lat)')} is used to indicate the response 'y' follows a spatial linear mixed model with a mean structure consisting of a linear combination of the two covariates 'X1' and 'X2' and a spatial random effect that has a Mat\'{e}rn covariance structure where the location coordinates are represented by 'lon' and 'lat' columns in \code{y_loc}. Note that when \code{gc_dist = TRUE}, the spatial random effect term should always take the form of 'Matern(1| lon + lat)'.
#' As discussed in the manuscript, it is assumed that the Mat\'{e}rn smoothness parameter for the spatial random effects is known, which can be supplied as \code{nu_rho}. The smoothness parameter can be chosen based on exploratory analysis of the data.
#' The argument \code{smooth_list} allows user to obtain the estimated conditional smoother if required. If \code{is.null(smooth_list)}, then no conditional smoother is estimated, otherwise the function will also estimate the conditional smoothers. If \code{smooth_list} is a list, then the elements of the list should have names matching those in \code{x_name}. Each element of the list is used to obtain a conditional smoother for a particular covariate. More specifically, each element of the list is a data frame with \code{length(x_name)} columns where the column names should match \code{x_name} and the column sharing the same name as the list element consists of a sequence of values supplied by the user while each of the other columns is a vector of chosen constant value. 
#' 
#' @param l Number of nearest-neighbours used to construct spatially aligned dataset.
#' @param dist_mat_y_x A distance matrix with \code{nrow(y_loc)} number of rows (for the response locations) and \code{nrow(x_loc)} number of columns (for the misaligned covariate locations), where the distance measure should be great-circle distance if \code{gc_dist = TRUE} or Euclidean distance if \code{gc_dist = FALSE}.
#' @param data_x A data frame consisting of the covariates at the misaligned covariate locations, with the column names matching \code{x_name} and \code{loc_name}.
#' @param data_y A data frame consisting of the response at the response location, with the column name matching \code{y_name} and \code{loc_name}.
#' @param gc_dist Logical. If \code{TRUE} great-circle distance (in kilometers) is used in the Mat\'{e}rn covariance matrices, otherwise Euclidean distance is used. 
#' @param dist_vec A vector consisting of the lower triangular part of the distance matrix, which can be obtained through \code{dist_mat[lower.tri(dist_mat)]} where \code{dist_mat} is a distance matrix. The details of \code{dist_mat} is given under 'Details'. 
#' @param y_name A character string representing the name of the response.
#' @param x_name A character vector containing the names of the covariates.
#' @param loc_name A character vector containing the names of the coordinates. When \code{gc_dist = TRUE}, it should be equal to \code{c('lon','lat')}.
#' @param slmm_formula A model formula describing the spatial linear mixed model of the response. The details of this is given under 'Details'. 
#' @param nu_rho Mat\'{e}rn smoothness parameter for the spatial random effect that is assumed to be known.
#' @param smooth_list Can be either \code{NULL} or a list of \code{length(x_name)} elements with names matching those in \code{x_name}. The details of this is given under 'Details'.
#' @param discretized_num Desired length of the sequence used to construct columns of varying values for each covariate. This only needs to be specified when \code{smooth_list = 'mean'}.
#' @param smooth_CI_level Confidence level for constructing the confidence interval of the conditional smoothers with the default being 0.95, and is only used when \code{!is.null(smooth_list)}.
#' 
#' @return A list with the elements being input arguments of the function call and the following additional elements:
#' \item{nearest_loc_mat:}{A matrix where each column represents the indices of l nearest neighbours for each location in \code{y_loc}.}
#' \item{hat_x_df:}{A data frame where each row consists of the mean of the covariate values at the l nearest neighbors for each location in \code{y_loc}. }
#' \item{fitme_fit:}{A \code{spaMM::fitme} model object based on fitting the spatial linear mixed model on the spatially aligned dataset.}
#' \item{fitme_fit_time:}{Time taken for the \code{spaMM::fitme} fit on the spatially aligned dataset.}
#' \item{fitme_warnings:}{A list of warnings for events that may have occurred during the \code{spaMM::fitme} fit on the spatially aligned dataset.}
#' \item{hatbeta:}{Estimated mean regression coefficients from the \code{spaMM::fitme} fit on the spatially aligned dataset.}
#' \item{hatsigma2_rho:}{Estimated variance parameter of the spatial random effect from the \code{spaMM::fitme} fit on the spatially aligned dataset.}
#' \item{hatnu_rho:}{Mat\'{e}rn smoothness parameter for the spatial random effect that is assumed known, which is the same as \code{nu_rho}.}
#' \item{hatalpha_rho:}{Estimated Mat\'{e}rn range parameter for the spatial random effect from the \code{spaMM::fitme} fit on the spatially aligned dataset.}
#' \item{hattau_epsilon:}{Estimated variance parameter of the residuals, or equivalently estimated nugget parameter for the response from the \code{spaMM::fitme} fit on the spatially aligned dataset.}
#' \item{hatSigma_err:}{Estimated covariance matrix of spatial random effect plus residuals.}
#' \item{naive_cov:}{Naive estimate of the covariance matrix for \code{hatbeta}, which ignore the uncertainty associated with the correction for spatial misalignment.}
#' \item{smoother:}{If \code{!is.null(smooth_list)}, then it is a list with \code{length(x_name)} number of elements, and each element of the list is a data frame consisting of \code{3 + length(x_name)} columns. In particular, these data frames are similar to the elements of \code{smooth_list}, with the difference being the additional first three columns in the elements of \code{smoother}. The additional columns correspond to the predicted mean component of the response excluding the spatial random effect(first column), the corresponding lower (second column) and upper (third column) confidence limits, where the prediction is obtained based on \code{fitme_fit} and the covariate values given in the remaining columns. If \code{is.null(smooth_list)}, then it is \code{NULL}.}

lnmr = function(l, dist_mat_y_x, data_x, data_y, gc_dist, dist_vec, y_name, x_name, loc_name, slmm_formula, nu_rho, smooth_list = NULL, discretized_num = NULL, smooth_CI_level = 0.95){
  m = nrow(data_x); n = nrow(data_y); K = length(x_name)
  
  #Finding l nearest neighbors for each location in y_loc
  nearest_loc_mat = matrix(apply(dist_mat_y_x, 1, FUN = function(d_vec){
    order(d_vec)[1:l]
  }), nrow = l, ncol = n, byrow = FALSE)
  
  #Compute the average covariate values based on the l nearest neighbors
  hat_x_df =  data.frame(t(apply(nearest_loc_mat,2, FUN = function(locs){
    if(l > 1){
      return(apply(data_x[locs,x_name], 2, mean))
    }else{
      return(as.matrix(data_x[locs,x_name])) 
    }
  })) )
  hat_x_df = data.frame(hat_x_df, data_y[, loc_name])
  colnames(hat_x_df) = c(x_name, loc_name)
  
  #Fit spatial linear mixed model
  if(gc_dist){
    fitme_fit = fitme(
      slmm_formula, 
      data = data.frame(data_y, hat_x_df[,x_name]),
      method="ML",
      control.dist=list(dist.method="Earth"),
      control = list(refit = TRUE),
      fixed = list(nu = nu_rho)
    )
  }else{
    fitme_fit = fitme(
      slmm_formula, 
      data = data.frame(data_y, hat_x_df[,x_name]),
      method="ML",
      control = list(refit = TRUE),
      fixed = list(nu = nu_rho)
    )
  }
  
  #Compute naive covariance matrix of beta coef
  hatSigma_err = as.numeric(fitme_fit$lambda) * vec2symMat(geoR::matern(u = dist_vec , kappa = get_ranPars(fitme_fit,which = 'corrPars')$'1'$nu, phi = (1/get_ranPars(fitme_fit,which = 'corrPars')$'1'$rho )), diag = F ) + fitme_fit$phi * diag(m+n)
  hatSigma_err = hatSigma_err[(m+1):(m+n),(m+1):(m+n)]
  naive_cov = solve( t(fitme_fit$X.pv) %*% solve(hatSigma_err) %*% fitme_fit$X.pv )
  
  #Compute smoother's and smoother's CI based on naive covariance matrix
  if(!is.null(smooth_list)){
    new_x_list = smooth_list
    
    if(!is.null(smooth_CI_level)){
      smoother = lapply(new_x_list, FUN = function(x_list){
        CI = get_intervals(object = fitme_fit, newdata = x_list, intervals = 'fixefVar', re.form  = NA, level = smooth_CI_level)
        df_out = data.frame(predict(object = fitme_fit, newdata = x_list, re.form = NA), CI , x_list)
        colnames(df_out) = c(y_name, paste0(y_name,'_lower'), paste0(y_name,'_upper'), x_name)
        df_out
      })
    }else{
      smoother = lapply(new_x_list, FUN = function(x_list){
        df_out = data.frame(predict(object = fitme_fit, newdata = x_list, re.form = NA), NA,NA ,x_list)
        colnames(df_out) = c(y_name, paste0(y_name,'_lower'), paste0(y_name,'_upper'), x_name)
        df_out
      })
    }
  }else{
    smoother = NULL
  }
  
  return(list(
    #Input Argument
    data_x = data_x, data_y = data_y, gc_dist = gc_dist, dist_mat_y_x = dist_mat_y_x, y_name = y_name, x_name = x_name, loc_name = loc_name, slmm_formula = slmm_formula, nu_rho = nu_rho, smooth_list = smooth_list, discretized_num = discretized_num, smooth_CI_level = smooth_CI_level,
    #Results
    nearest_loc_mat = nearest_loc_mat, hat_x_df = hat_x_df,
    fitme_fit = fitme_fit, fitme_fit_time = fitme_fit$how$fit_time, fitme_warnings = fitme_fit$warnings,
    hatbeta = fitme_fit$fixef, hatsigma2_rho = fitme_fit$lambda, hatnu_rho = get_ranPars(fitme_fit,which = 'corrPars')$'1'$nu, hatalpha_rho = get_ranPars(fitme_fit,which = 'corrPars')$'1'$rho, hattau_epsilon = fitme_fit$phi, 
    hatSigma_err = hatSigma_err, naive_cov = naive_cov, smoother = smoother
  ))
}