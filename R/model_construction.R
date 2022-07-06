#' Compute the precision matrix of the RW2 part:
#' 
#' @param x A vector of covariate values.
#' @param r An integer number of knots.
#' @param diag_noise A positive number denotes the diagonal adjustment.
#' @return A sparse matrix denotes the precision matrix of the RW2 component.
RW2_prec_comp <- function(x, r, diag_noise){
  d <- (r)
  a <- min(x)
  b <- max(x)
  b_basis <- fda::create.bspline.basis(rangeval = c(a,b), nbasis = r, norder = 4)
  coefs <- diag(rep(1, r))
  bfd <- fda::fd(coef = coefs, basisobj = b_basis)
  as(fda::inprod.bspline(fdobj1 = bfd, fdobj2 = bfd, nderiv1 = 2, nderiv2 = 2) + diag(diag_noise, ncol = r, nrow = r),"dgTMatrix")
}


#' Compute the design matrix of the RW2 part:
#' 
#' @param x A vector of covariate values.
#' @param r An integer number of knots.
#' @return A sparse matrix denotes the design matrix of the RW2 component.
RW2_design_comp <- function(x, r){
  d <- (r)
  a <- min(x)
  b <- max(x)
  b_basis <- fda::create.bspline.basis(rangeval = c(a,b), nbasis = r, norder = 4)
  as(fda::getbasismatrix(evalarg = x, basisobj = b_basis), "dgTMatrix")
}



#' Do the Approximate Bayesian Inference for the specified Cox 
#' Proportional Hazard Model.
#' 
#' @param data A dataframe that should contain all the variables listed as times, cens, fixed,
#' frailty and RW2.
#' @param times A single string of variable name to be used as survival times.
#' @param cens A single string of variable name to be used as censoring indicators
#' @param fixed A vector of strings of covariate names to be used as linear
#' fixed effects in the model. All covariates are assumed to be numeric.
#' @param frailty A single string of covariate name to be used as the group for
#' frailty in the model, default is NULL meaning no frailty.
#' @param RW2 A single string of covariate name to be used as the RW2 smoothing
#' covariate, default is NULL meaning no RW2 smoothing. The covariate is assumed to be
#' numeric.
#' @param fixed_control A list that specifies the prior precision of all the fixed
#' effect parameters. Default being list(betaprec = .001).
#' @param frailty_control A list that specifies the PC prior for the frailty SD
#' parameters. Default being list(alpha = 0.5, u = 1)
#' @param RW2_control A list that specifies the PC prior for the RW2 SD
#' parameters and number of knots used. Default being list(alpha = 0.5, u = 1, r = 50)
#' @param Inference_control A list that specifies the setting of the AGHQ inference, including
#' the number of grid points K. Default being list(aghq_k = 4).
#' @param diag_noise The small amount of diagonal adjustment to add when RW2 exists in the model.
#' The default value is diag_noise = 0.0001.
#' @return A list that contains the fitted AGHQ object of the model and other components.
#'  
#' @export 
abcoxp_fit <- function(data, times, cens, fixed = NULL, frailty = NULL, RW2 = NULL, 
                   fixed_control = list(betaprec = .001), frailty_control = list(alpha = 0.5, u = 1), 
                   RW2_control = list(alpha = 0.5, u = 1, r = 50), 
                   Inference_control = list(aghq_k = 4), diag_noise = 0.0001){
  data <- cbind(data[times], data[cens], data[fixed], data[frailty], data[RW2])
  colnames(data) <- c("times", "cens", fixed, frailty, RW2)
  data <- dplyr::arrange(data, by = times)
  data$ranks <- rank(data$times, ties.method = "min")
  n <- nrow(data)
  D <- cbind(Matrix::Matrix(1,n-1,1),Matrix::Diagonal(n-1,-1))
  if(length(fixed) > 0 & length(frailty) == 1 & length(RW2) == 1){
    ## Create X, B1, B2 design
    ## Create P2 precision
    X <- as(as.matrix(data[fixed]), "dgTMatrix")
    u1 <- as.numeric(data[frailty][,1])
    B1 <- Matrix::Diagonal(n = n)[match(u1,unique(u1)),order(unique(u1))]
    u2 <- as.numeric(data[RW2][,1])
    B2 <- RW2_design_comp(x = u2, r = RW2_control$r)
    P2 <- RW2_prec_comp(x = u2, r = RW2_control$r, diag_noise = diag_noise)
    tmbdat <- list(
      # Design matrix
      X = X,
      B1 = B1,
      B2 = B2,
      # Penalty matrix
      P2 = P2,
      # Differencing matrix
      D = D,
      # Log determinant of penalty matrix (without the sigma part)
      logP2det = as.numeric(Matrix::determinant(P2,logarithm = TRUE)$modulus),
      # Response
      ranks = as.integer(data$ranks),
      cens = as.integer(data$cens),
      # Prior params
      u1 = frailty_control$u,
      alpha1 = frailty_control$alpha,
      u2 = RW2_control$u,
      alpha2 = RW2_control$alpha,
      betaprec = fixed_control$betaprec
    )
    tmbparams <- list(
      W = rep(0,ncol(B1) + ncol(B2) + ncol(X)),
      theta1 = 0, # -2log(sigma)
      theta2 = 0
    )
    ff <- TMB::MakeADFun(
      data = c(model = "fixed_smooth_frailty", tmbdat),
      parameters = tmbparams,
      random = "W",
      DLL = "abcoxp_TMBExports",
      silent = TRUE
    )
    # Hessian not implemented for RE models
    ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    # AGHQ
    quad <- aghq::marginal_laplace_tmb(ff, k = Inference_control$aghq_k, c(0,0))
    return(list(model = quad, data = data, X = X, B1 = B1, B2 = B2, components = list(fixed = fixed, frailty = frailty, RW2 = RW2)))
  }
  else if(length(fixed) > 0 & length(frailty) == 1 & length(RW2) == 0){
    ## Create X, B1 design
    X <- as(as.matrix(data[fixed]), "dgTMatrix")
    u1 <- as.numeric(data[frailty][,1])
    B1 <- Matrix::Diagonal(n = n)[match(u1,unique(u1)),order(unique(u1))]
    tmbdat <- list(
      # Design matrix
      X = X,
      B1 = B1,
      # Differencing matrix
      D = D,
      # Response
      ranks = as.integer(data$ranks),
      cens = as.integer(data$cens),
      # Prior params
      u = frailty_control$u,
      alpha = frailty_control$alpha,
      betaprec = fixed_control$betaprec
    )
    tmbparams <- list(
      W = rep(0,ncol(B1) + ncol(X)),
      theta = 0 # -2log(sigma)
    )
    ff <- TMB::MakeADFun(
      data = c(model = "fixed_frailty", tmbdat),
      parameters = tmbparams,
      random = "W",
      DLL = "abcoxp_TMBExports",
      silent = TRUE
    )
    # Hessian not implemented for RE models
    ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    # AGHQ
    quad <- aghq::marginal_laplace_tmb(ff, k = Inference_control$aghq_k, c(0))
    return(list(model = quad, data = data, X = X, B1 = B1, components = list(fixed = fixed, frailty = frailty, RW2 = RW2)))
  }
  else if(length(fixed) > 0 & length(frailty) == 0 & length(RW2) == 1){
    ## Create X, B2 design
    ## Create P2 precision
    X <- as(as.matrix(data[fixed]), "dgTMatrix")
    u2 <- as.numeric(data[RW2][,1])
    B2 <- RW2_design_comp(x = u2, r = RW2_control$r)
    P2 <- RW2_prec_comp(x = u2, r = RW2_control$r, diag_noise = diag_noise)
    tmbdat <- list(
      # Design matrix
      X = X,
      B2 = B2,
      # Penalty matrix
      P2 = P2,
      # Differencing matrix
      D = D,
      # Log determinant of penalty matrix (without the sigma part)
      logP2det = as.numeric(Matrix::determinant(P2,logarithm = TRUE)$modulus),
      # Response
      ranks = as.integer(data$ranks),
      cens = as.integer(data$cens),
      # Prior params
      u = RW2_control$u,
      alpha = RW2_control$alpha,
      betaprec = fixed_control$betaprec
    )
    tmbparams <- list(
      W = rep(0,ncol(B2) + ncol(X)),
      theta = 0
    )
    ff <- TMB::MakeADFun(
      data = c(model = "fixed_smooth", tmbdat),
      parameters = tmbparams,
      random = "W",
      DLL = "abcoxp_TMBExports",
      silent = TRUE
    )
    # Hessian not implemented for RE models
    ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    # AGHQ
    ## Here with bug currently
    quad <- aghq::marginal_laplace_tmb(ff, k = Inference_control$aghq_k, startingvalue = 0)
    return(list(model = quad, data = data, X = X, B2 = B2, components = list(fixed = fixed, frailty = frailty, RW2 = RW2)))
  }
  else if(length(fixed) > 0 & length(frailty) == 0 & length(RW2) == 0){
    ## Create X design
    X <- as(as.matrix(data[fixed]), "dgTMatrix")
    tmbdat <- list(
      # Design matrix
      X = X,
      # Differencing matrix
      D = D,
      # Response
      ranks = as.integer(data$ranks),
      cens = as.integer(data$cens),
      # Prior params
      betaprec = fixed_control$betaprec
    )
    tmbparams <- list(
      W = rep(0,ncol(X))
    )
    ff <- TMB::MakeADFun(
      data = c(model = "fixed", tmbdat),
      parameters = tmbparams,
      # random = "W",
      DLL = "abcoxp_TMBExports",
      silent = T
    )
    # Hessian not implemented for RE models
    ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    opt <- nlminb(start = ff$par, objective = ff$fn, gradient = ff$gr, hessian = ff$he, 
                  control = list(eval.max = 20000, iter.max = 20000))
    prec_matrix <- Matrix::forceSymmetric(ff$he(opt$par))
    return(list(model = list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt), data = data, X = X, components = list(fixed = fixed, frailty = frailty, RW2 = RW2)))
  }
  else if(length(fixed) == 0 & length(frailty) == 0 & length(RW2) == 1){
    ## Create B2 design
    ## Create P2 precision
    u2 <- as.numeric(data[RW2][,1])
    B2 <- RW2_design_comp(x = u2, r = RW2_control$r)
    P2 <- RW2_prec_comp(x = u2, r = RW2_control$r, diag_noise = diag_noise)
    tmbdat <- list(
      # Design matrix
      B2 = B2,
      # Penalty matrix
      P2 = P2,
      # Differencing matrix
      D = D,
      # Log determinant of penalty matrix (without the sigma part)
      logP2det = as.numeric(Matrix::determinant(P2,logarithm = TRUE)$modulus),
      # Response
      ranks = as.integer(data$ranks),
      cens = as.integer(data$cens),
      # Prior params
      u = RW2_control$u,
      alpha = RW2_control$alpha
    )
    tmbparams <- list(
      W = rep(0, ncol(B2)),
      theta = 0
    )
    ff <- TMB::MakeADFun(
      data = c(model = "smooth", tmbdat),
      parameters = tmbparams,
      random = "W",
      DLL = "abcoxp_TMBExports",
      silent = TRUE
    )
    # Hessian not implemented for RE models
    ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    # AGHQ
    quad <- aghq::marginal_laplace_tmb(ff, k = Inference_control$aghq_k, c(0))
    return(list(model = quad, data = data, B2 = B2, components = list(fixed = fixed, frailty = frailty, RW2 = RW2)))
  }
  else{
    stop("Model does not have a TMB template yet, need to manually implement.")
  }
  
}



