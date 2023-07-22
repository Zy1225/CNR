#' @title 
#' Cokrig-and-regress for spatially misaligned data
#' 
#' @description 
#' Perform cokrig-and-regress on spatially misaligned data.
#' 
#' @details 
#' This function perform cokrig-and-regress (CNR) on spatially misaligned data, where the continuous response is spatially misaligned with multiple correlated continuous covariates. 
#' As discussed in the manuscript, it is assumed that the Mat\'{e}rn smoothness parameters for the covariates and the spatial random effects are known, which can be supplied as \code{nu_vec_x} and \code{nu_rho}. These smoothness parameters can be chosen based on exploratory analysis of the data.
#' The distance matrix \code{dist_mat} should be constructed by first combining \code{x_loc} and \code{y_loc} (with \code{x_loc} ordered first) and then computing the distance matrix based on all locations, where the distance measure should be great-circle distance if \code{gc_dist = TRUE} or Euclidean distance if \code{gc_dist = FALSE}.
#' The model formula \code{slmm_formula} should be a linear model formula that can be handled by \code{spaMM::fitme} fitting function with a spatial random effect term that is assumed to have a Mat\'{e}rn covariance structure. For instance, \code{slmm_formula = as.formula('y ~ X1 + X2 + Matern(1| lon + lat)')} is used to indicate the response 'y' follows a spatial linear mixed model with a mean structure consisting of a linear combination of the two covariates 'X1' and 'X2' and a spatial random effect that has a Mat\'{e}rn covariance structure where the location coordinates are represented by 'lon' and 'lat' columns in \code{y_loc}. Note that when \code{gc_dist = TRUE}, the spatial random effect term should always take the form of 'Matern(1| lon + lat)'.
#' 
#'  The argument \code{smooth_list} allows user to obtain the estimated conditional smoother if required. If \code{is.null(smooth_list)}, then no conditional smoother is estimated, otherwise the function will also estimate the conditional smoothers. If \code{smooth_list} is a list, then the elements of the list should have names matching those in \code{x_name}. Each element of the list is used to obtain a conditional smoother for a particular covariate. More specifically, each element of the list is a data frame with \code{length(x_name)} columns where the column names should match \code{x_name} and the column sharing the same name as the list element consists of a sequence of values supplied by the user while each of the other columns is a vector of chosen constant value. If \code{smooth_list = 'mean}, then a similar list is constructed where for the list element with name \code{x_name[k]}, the sequence of values of the column with name \code{x_name[k]} is set to be a regular sequence of \code{discretized_num} number of values from the smallest to the largest values of \code{data_x[, x_name[k]]} while the other columns are set to be their estimated marginal mean parameters i.e., values from \code{hatmu_vec_x} in 'Value'.
#' 
#' @param data_x A data frame consisting of the covariates at the misaligned covariate locations, with the column names matching \code{x_name} and \code{loc_name}.
#' @param data_y A data frame consisting of the response at the response location, with the column name matching \code{y_name} and \code{loc_name}.
#' @param gc_dist Logical. If \code{TRUE} great-circle distance (in kilometers) is used in the Mat\'{e}rn covariance matrices, otherwise Euclidean distance is used. 
#' @param nu_vec_x A vector of Mat\'{e}rn smoothness parameters for the covariates, assumed to be known.
#' @param dist_vec A vector consisting of the lower triangular part of the distance matrix, which can be obtained through \code{dist_mat[lower.tri(dist_mat)]} where \code{dist_mat} is a distance matrix. The details of \code{dist_mat} is given under 'Details'. 
#' @param y_name A character string representing the name of the response.
#' @param x_name A character vector containing the names of the covariates.
#' @param loc_name A character vector containing the names of the coordinates. When \code{gc_dist = TRUE}, it should be equal to \code{c('lon','lat')}.
#' @param slmm_formula A model formula describing the spatial linear mixed model of the response. The details of this is given under 'Details'. 
#' @param nu_rho Mat\'{e}rn smoothness parameter for the spatial random effect that is assumed to be known. The details of this is given under 'Details'.
#' @param smooth_list Can be one of \code{NULL}, a character vector \code{'mean}, or a list of \code{length(x_name)} elements with names matching those in \code{x_name}. 
#' @param discretized_num Desired length of the sequence used to construct columns of varying values for each covariate. This only needs to be specified when \code{smooth_list = 'mean'}.
#' @param smooth_CI_level Confidence level for constructing the confidence interval of the conditional smoothers with the default being 0.95, and is only used when \code{!is.null(smooth_list)}.
#' 
#' @return A list with the elements being input arguments of the function call and the following additional elements:
#' \item{step1_warnings:}{A list of \code{length(x_name)} number of elements where each element is a list of warnings for events that may have occurred during the \code{spaMM::fitme} fit in Step 1 of CNR. }
#' \item{step1_fitme_time_total:}{Total time taken for all \code{spaMM::fitme} fits in Step 1 of CNR. }
#' \item{hatnu_vec_x:}{A vector of Mat\'{e}rn smoothness parameters for the covariates that are assumed known, which is the same as \code{nu_vec_x} in 'Arguments'.}
#' \item{hatalpha_vec_x:}{A vector of estimated Mat\'{e}rn range parameters for the covariates from Step 1 of CNR.}
#' \item{hatsigma2_vec_x:}{A vector of estimated Mat\'{e}rn variance parameters for the covariates from Step 1 of CNR.}
#' \item{hattau_vec_x:}{A vector of estimated nugget parameters for the covariates from Step 1 of CNR.}
#' \item{hatmu_vec_x:}{A vector of estimated mean parameters for the covariates from Step 1 of CNR.}
#' \item{hatR:}{Estimated unstructured cross-correlation matrix from Step 1 of CNR.}
#' \item{hat_x_df:}{A data frame consisting of the cokriging prediction of the covariates from Step 2 of CNR.}
#' \item{step3_fitme:}{A \code{spaMM::fitme} model object based on fitting the spatial linear mixed model in Step 3 of CNR.}
#' \item{step3_fitme_time:}{Time taken for the \code{spaMM::fitme} fit in Step 3 of CNR.}
#' \item{step3_warnings:}{A list of warnings for events that may have occurred during the \code{spaMM::fitme} fit in Step 3 of CNR.}
#' \item{hatbeta:}{Estimated mean regression coefficients from Step 3 of CNR.}
#' \item{hatsigma2_rho:}{Estimated variance parameter of the spatial random effect from Step 3 of CNR.}
#' \item{hatnu_rho:}{Mat\'{e}rn smoothness parameter for the spatial random effect that is assumed known, which is the same as \code{nu_rho}.}
#' \item{hatalpha_rho:}{Estimated Mat\'{e}rn range parameter for the spatial random effect from Step 3 of CNR.}
#' \item{hattau_epsilon:}{Estimated variance parameter of the residuals, or equivalently estimated nugget parameter for the response.}
#' \item{hatSigma_err:}{Estimated covariance matrix of spatial random effect plus residuals.}
#' \item{naive_cov:}{Naive estimate of the covariance matrix for \code{hatbeta}, which ignore the uncertainty associated with the cokriging prediction of covariates.}
#' \item{smooth_list_out:}{Equals to \code{smooth_list} if \code{smooth_list} if a list, or the constructed list described in 'Details' if \code{smooth_list = 'mean'}. Otherwise, equals to \code{NULL} if \code{is.null(smooth_list)}.}
#' \item{smoother:}{If \code{!is.null(smooth_list)}, then it is a list with \code{length(x_name)} number of elements, and each element of the list is a data frame consisting of \code{3 + length(x_name)} columns. In particular, these data frames are similar to the elements of \code{smooth_list} (or the constructed list when \code{smooth_list = 'mean'}), with the difference being the additional first three columns in the elements of \code{smoother}. The additional colums correspond to the predicted mean component of the response excluding the spatial random effect  (first column), the corresponding lower (second column) and upper (third column) confidence limits, where the prediction is obtained based on \code{step3_fitme} and the covariate values given in the remaining columns. If \code{is.null(smooth_list)}, then it is \code{NULL}.}


cnr = function(data_x, data_y, x_loc, y_loc, gc_dist, nu_vec_x, dist_vec, y_name, x_name, loc_name, slmm_formula, nu_rho, smooth_list = NULL, discretized_num = NULL, smooth_CI_level = 0.95){
  m = nrow(data_x); n = nrow(data_y); K = length(x_name)
  
  #CNR Step 1
  if(gc_dist){
    step1_fitme_list = lapply(1:K, function(k){
      fitme(formula = as.formula(paste(x_name[k] , paste('1 + Matern(1|lon+lat)') , sep = " ~ ")),
            data = data.frame(data_x),
            method = 'ML',
            control.dist=list(dist.method="Earth"),
            control = list(refit = TRUE),
            fixed = list(nu = nu_vec_x[k])
      )
    })
  }else{
    step1_fitme_list = lapply(1:K, function(k){
      fitme(formula = as.formula(paste(x_name[k] , paste('1 + Matern(1|', paste(loc_name, collapse = '+') ,')') , sep = " ~ ")),
            data = data.frame(data_x),
            method = 'ML',
            control = list(refit = TRUE),
            fixed = list(nu = nu_vec_x[k])
      )
    })
  }
  step1_warnings = sapply(step1_fitme_list,function(x){x$warnings})
  step1_fitme_time_total = sum(sapply(step1_fitme_list, FUN = function(x){x$how$fit_time}))
  
  hatnu_vec_x = sapply(step1_fitme_list, FUN = function(x){get_ranPars(x,which = 'corrPars')$'1'$nu})
  hatalpha_vec_x = sapply(step1_fitme_list, FUN = function(x){get_ranPars(x,which = 'corrPars')$'1'$rho})
  hatsigma2_vec_x = sapply(step1_fitme_list, FUN = function(x){x$lambda})
  hattau_vec_x = sapply(step1_fitme_list, FUN = function(x){x$phi})
  hatmu_vec_x = sapply(step1_fitme_list, FUN = function(x){x$fixef})
  
  hatR = (1/m) * crossprod(sapply(1:K, FUN = function(k){
    hatSigma_k = hatsigma2_vec_x[k]*vec2symMat(geoR::matern(u = dist_vec, kappa = hatnu_vec_x[k], phi = (1/(hatalpha_vec_x[k]))), diag = F) + hattau_vec_x[k] * diag(m+n)
    hatL_k_tildeS_T = chol(hatSigma_k[1:m, 1:m])
    crossprod( solve(hatL_k_tildeS_T), as.matrix(as.matrix(data_x[,x_name[k]]) - rep(hatmu_vec_x[k],m)) )
    # crossprod( solve(hatL_k_tildeS_T), data_x[,x_name[k]] - rep(hatmu_vec_x[k],m))
  }) )
  
  #CNR Step 2
  hatSigma_k_list = lapply(1:K, function(k){
    as.numeric(hatsigma2_vec_x[k]) * vec2symMat(geoR::matern(u = dist_vec , kappa = hatnu_vec_x[k], phi = (1/(hatalpha_vec_x[k] ))), diag = F ) + hattau_vec_x[k] * diag(m+n)
  })
  
  hat_x_df = data.frame(sapply(1:K, function(k){
    rep(hatmu_vec_x[k], n) + hatSigma_k_list[[k]][ (m+1):(m+n)  , 1:m ] %*% solve(hatSigma_k_list[[k]][ 1:m, 1:m ]) %*% ( data_x[,x_name[k]] - rep(hatmu_vec_x[k],m)  )
  })  )
  hat_x_df = data.frame(hat_x_df, data_y[, loc_name])
  colnames(hat_x_df) = c(x_name, loc_name)
  
  #CNR Step 3
  if(gc_dist){
    step3_fitme = fitme(
      slmm_formula, 
      data = data.frame(data_y, hat_x_df[,x_name]),
      method="ML",
      control.dist=list(dist.method="Earth"),
      control = list(refit = TRUE),
      fixed = list(nu = nu_rho)
    )
  }else{
    step3_fitme = fitme(
      slmm_formula, 
      data = data.frame(data_y, hat_x_df[,x_name]),
      method="ML",
      control = list(refit = TRUE),
      fixed = list(nu = nu_rho)
    )
  }
  
  #Compute naive covariance matrix of beta coefficient estimates
  hatSigma_err = as.numeric(step3_fitme$lambda) * vec2symMat(geoR::matern(u = dist_vec , kappa = get_ranPars(step3_fitme,which = 'corrPars')$'1'$nu, phi = (1/get_ranPars(step3_fitme,which = 'corrPars')$'1'$rho )), diag = F ) + step3_fitme$phi * diag(m+n)
  hatSigma_err = hatSigma_err[(m+1):(m+n),(m+1):(m+n)]
  naive_cov = solve( t(step3_fitme$X.pv) %*% solve(hatSigma_err) %*% step3_fitme$X.pv )
  
  
  #Compute smoother's and smoother's CI based on naive covariance matrix
  if(!is.null(smooth_list)){
    if(smooth_list[1] == 'mean'){
      new_x_list = lapply(1:K, FUN = function(k){
        temp_df = matrix(0, nrow = discretized_num, ncol = K, dimnames = list(c(),x_name))
        temp_df[,x_name[k]] = seq(from = min(data_x[,x_name[k]]), to = max(data_x[,x_name[k]]), length.out = discretized_num)
        temp_df[, x_name[-k] ] = rep(hatmu_vec_x[-k], each = discretized_num)
        temp_df
      })
      names(new_x_list) = x_name
      smooth_list_out = new_x_list
    }else{
      new_x_list = smooth_list
      smooth_list_out = smooth_list
    }
    
    if(!is.null(smooth_CI_level)){
      smoother = lapply(1:K, FUN = function(k){
        CI = get_intervals(object = step3_fitme, newdata = new_x_list[[x_name[k]]], intervals = 'fixefVar', re.form  = NA, level = smooth_CI_level)
        df_out = data.frame(predict(object = step3_fitme, newdata = new_x_list[[x_name[k]]], re.form = NA), CI , new_x_list[[x_name[k]]])
        colnames(df_out)[1:3] = c(y_name, paste0(y_name,'_lower'), paste0(y_name,'_upper'))
        df_out
      })
    }else{
      smoother = lapply(1:K, FUN = function(k){
        df_out = data.frame(predict(object = step3_fitme, newdata = new_x_list[[x_name[k]]], re.form = NA), NA,NA ,new_x_list[[x_name[k]]])
        colnames(df_out)[1:3] = c(y_name, paste0(y_name,'_lower'), paste0(y_name,'_upper'))
        df_out
      })
    }
  }else{
    smoother = NULL
    smooth_list_out = NULL
  }
  
  return(list(
    #Input Argument
    data_x = data_x, data_y = data_y, gc_dist = gc_dist, nu_vec_x = nu_vec_x, dist_vec = dist_vec, y_name = y_name, x_name = x_name, loc_name = loc_name, slmm_formula = slmm_formula, nu_rho = nu_rho, smooth_list = smooth_list, discretized_num = discretized_num, smooth_CI_level = smooth_CI_level,
    #Step 1 results
    step1_warnings = step1_warnings, step1_fitme_time_total = step1_fitme_time_total, hatnu_vec_x = hatnu_vec_x, hatalpha_vec_x = hatalpha_vec_x, hatsigma2_vec_x = hatsigma2_vec_x, hattau_vec_x = hattau_vec_x, hatmu_vec_x = hatmu_vec_x, hatR = hatR,
    #Step 2 results
    hat_x_df = hat_x_df,
    #Step 3 results
    step3_fitme = step3_fitme, step3_fitme_time = step3_fitme$how$fit_time, step3_warnings = step3_fitme$warnings,
    hatbeta = step3_fitme$fixef, hatsigma2_rho = step3_fitme$lambda, hatnu_rho = get_ranPars(step3_fitme,which = 'corrPars')$'1'$nu, hatalpha_rho = get_ranPars(step3_fitme,which = 'corrPars')$'1'$rho, hattau_epsilon = step3_fitme$phi, 
    hatSigma_err = hatSigma_err, naive_cov = naive_cov, smooth_list_out = smooth_list_out, smoother = smoother
  ))
}
