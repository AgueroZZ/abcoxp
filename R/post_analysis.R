#' Implements the Approximate Bayesian Inference for the specified Cox 
#' Proportional Hazard Model with Analysis of Posterior.
#' 
#' @param fitted_model A list from the output of `abcoxp_fit`.
#' @param n_samps An integer value of the number of samples.
#' @return A list that contains the posterior samples.
#' 
abcoxp_sampling <- function(fitted_model, nsamps){
  if(length(fitted_model$components$RW2) == 0 & length(fitted_model$components$frailty) == 0){
    beta_samps <- LaplacesDemon::rmvnp(n = nsamps, mu = fitted_model$model$mean, Omega = as.matrix(fitted_model$model$prec))
    colnames(beta_samps) <- c(fitted_model$components$fixed)
    return(list(fixed_samps = beta_samps))
  }
  else{
    all_samps <- aghq::sample_marginal(fitted_model$model, M = nsamps)
    if(length(fitted_model$components$RW2)== 1 & length(fitted_model$components$frailty) == 1){
      B1 <- fitted_model$B1
      B2 <- fitted_model$B2
      K1 <- ncol(B1)
      K2 <- ncol(B2)
      frailty_samps <- t(all_samps$samps[1:K1,,drop=FALSE])
      RW2_weights_samps <- t(all_samps$samps[(K1+1):(K1+K2),,drop=FALSE])
      RW2_values_samps <- cbind(x = fitted_model$data[[fitted_model$components$RW2]], B2 %*% all_samps$samps[(K1+1):(K1+K2),])
      RW2_values_samps <- as.data.frame(as.matrix(RW2_values_samps))
      RW2_values_samps <- dplyr::distinct(dplyr::arrange(RW2_values_samps, by = x), x, 
                                          .keep_all = TRUE)
      beta_samps <- t(all_samps$samps[-c(1:(K1+K2)),,drop=FALSE])
      colnames(beta_samps) <- c(fitted_model$components$fixed)
      return(list(fixed_samps = beta_samps, frailty_samps = frailty_samps, 
                  RW2_values_samps = RW2_values_samps, RW2_weights_samps = RW2_weights_samps,
                  frailty_SD = sqrt(1/exp(all_samps$theta1)),
                  RW2_SD = sqrt(1/exp(all_samps$theta2))))
      
    }
    else if(length(fitted_model$components$RW2)== 0 & length(fitted_model$components$frailty) == 1){
      B1 <- fitted_model$B1
      K1 <- ncol(B1)
      K2 <- 0
      frailty_samps <- t(all_samps$samps[1:K1,,drop=FALSE])
      beta_samps <- t(all_samps$samps[-c(1:(K1+K2)),,drop=FALSE])
      colnames(beta_samps) <- c(fitted_model$components$fixed)
      return(list(fixed_samps = beta_samps, frailty_samps = frailty_samps, 
                  frailty_SD = sqrt(1/exp(all_samps$theta))))
      
    }
    else if(length(fitted_model$components$RW2)== 1 & length(fitted_model$components$frailty) == 0){
      B2 <- fitted_model$B2
      K2 <- ncol(B2)
      K1 <- 0
      RW2_weights_samps <- t(all_samps$samps[(K1+1):(K1+K2),,drop=FALSE])
      RW2_values_samps <- cbind(x = fitted_model$data[[fitted_model$components$RW2]], B2 %*% all_samps$samps[(K1+1):(K1+K2),])
      RW2_values_samps <- as.data.frame(as.matrix(RW2_values_samps))
      RW2_values_samps <- dplyr::distinct(dplyr::arrange(RW2_values_samps, by = x), x, 
                                          .keep_all = TRUE)
      beta_samps <- t(all_samps$samps[-c(1:(K1+K2)),,drop=FALSE])
      colnames(beta_samps) <- c(fitted_model$components$fixed)
      return(list(fixed_samps = beta_samps, 
                  RW2_values_samps = RW2_values_samps, RW2_weights_samps = RW2_weights_samps,
                  RW2_SD = sqrt(1/exp(all_samps$theta))))
      
    }
  }
}






