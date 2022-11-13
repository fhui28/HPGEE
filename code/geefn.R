#' @title Fit unpenalized GEEs
#' 
#' @description Function specifically used to fit unpenalized GEEs, assuming either an independence, an unstructured, or a reduced-rank working correlation matrix.
#' 
#' @param y A response matrix i.e., multivariate abundance data.
#' @param X Model matrix. It is *always* assumed that the first column in \code{X} is a vector of ones corresponding to an intercept term. 
#' @param bigX An optional enlarged model matrix that is based on \code{X}. This can be passed in if \code{bigX} has was created externally for other reasons. Otherwise, it is constructed internally within the function.
#' @param family The family to use, which then implies the link function and mean-variance relationship to use in the GEE.
#' @param corstr The structure to use for the working correlation in the GEE. Currently, the options permitted are "independence", "unstructured", and "reducedrank".
#' @param rank The rank of the reduced-rank working correlation matrix, if \code{corstr = "reducedrank"}.
#' @param start_coef A *vector* of starting values for the coefficients. Note this has be supplied as a vector (with regression coefficients running on a per-response basis) matrix. 
#' @param num_cores The number of cores to use as part of computation. This only comes into play during the initial set up of the quantities to pass into the GEE. Defaults to \code{doParallel::detectCores() - 1}.
#' @param method Method to use for updating regression coefficients. Current choices are "glmnet" for coordinate descent via the \code{glmnet} package (default), and "newtonraphson" for Newton-Raphson.
#' @param control A list for controlling the GEE estimation algorithm, including the maximum number of iterations (maxit), tolerance (tol), whether to employ step halving when updating the regression coefficients (step_halving), maximum number of step_halving attempts (max_halving_steps), and whether to print updates along the way (trace)
#' 
#' @return A object of class "unpenalized_gee" with the following values, as appropriate:
#' \item{call:}{The function call.}
#' \item{coefficients:}{The full vector of regression coefficients (running on a per-response basis).}
#' \item{coefficients_mattrix:}{The regression coefficients but arranged into a matrix, with the number of rows equal to the number of responses i.e., \code{ncol(y)}.}
#' \item{working_correlation:}{The working correlation matrix.}
#' \item{disp_param:}{The response-specific dispersion parameters. Note families which do not require dispersion parameters will still return a vector of values, although these are not updated/can be ignored.}
#' \item{power_param:}{The response-specific power parameters. Note families which do not require power parameters will still return a vector of values, although these are not updated/can be ignored.}
#' \item{linear_predictor:}{The matrix of linear predictors. The number of columns is equal to the number of responses i.e., \code{ncol(y)}.}
#' \item{fitted_values:}{The matrix of fitted values. The number of columns is equal to the number of responses i.e., \code{ncol(y)}.}
#' \item{residuals:}{The matrix of raw residuals. The number of columns is equal to the number of responses i.e., \code{ncol(y)}.}
#' \item{sandwich_matrix:}{The estimated sandwich covariance matrix for the full vector of regression coefficients i.e., corresponding to \code{coefficients} }
#' \item{stderr_matrix:}{An estimated matrix of standard errors for the regression coefficients i.e., corresponding to \code{coefficient_matrix}, with the number of rows equal to the number of responses i.e., \code{ncol(y)}.}


## TODO: Speed up estimation even more...most time taken in glmnet though?

library(doParallel)
library(foreach)
library(glmnet)
library(MASS)
library(Matrix)
#library(ncvreg)
#library(svMisc)
library(tweedie)

function() {
   y = simy
   family = myfam
   bigX = NULL
   control = list(trace = TRUE)
   corstr = "reducedrank"
   method = "glmnet"
   rank = 3
   num_cores = 1
   start_coef = NULL
   }

unpenalized_gee <- function(y, X, bigX = NULL, family = gaussian(), corstr = "independence", rank = 2, start_coef = NULL, num_cores = NULL, method = "glmnet",
                  control = list(maxit = 100, tol = 1e-6, step_halving = TRUE, max_halving_steps = 30, trace = FALSE)) {
   
   ##-----------------------
   ## Opening checks and set up
   ##-----------------------
   if(is.null(num_cores))
      num_cores <- detectCores() - 1
   registerDoParallel(cores = num_cores)
   
   method <- match.arg(method, choices = c("newtonraphson","glmnet")) 
   
   num_units <- nrow(y) 
   num_resp <- ncol(y) 
   num_X <- ncol(X)
   if(is.null(colnames(y)))
      colnames(y) <- paste0("resp", 1:num_resp)
   
   yvec <- as.vector(t(y)) # Response runs faster than unit
   # if(!is.null(offset))
   #    offset <- as.vector(t(offset))
   # if(is.null(offset))
   #    offset <- as(rep(0,length(yvec)),"sparseMatrix")
   X <- as.matrix(X)
   if(is.null(bigX)) {
      bigX <- NULL
      bigX <- foreach(i = 1:num_units, .combine = "rbind") %dopar% kronecker(Diagonal(n = num_resp),t(X[i,]))
      }

   corstr <- match.arg(corstr, choices = c("independence","unstructured", "reducedrank")) 
   if(corstr == "reducedrank") {
      if(num_resp < 3) 
         stop("A reduced rank correlation matrix requires at least three responses.")
    
      dof <- 0.5 * ((num_resp - rank)^2 - num_resp - rank)
      if(dof < 0)
         stop(sprintf(ngettext(rank,
                               "%d factor is too many for %d variables",
                               "%d factors are too many for %d variables"),
                      rank, num_resp), domain = NA)
      
      }
   if(!family$family %in% c("gaussian","binomial", "poisson", "negative_binomial", "tweedie", "Beta"))
      stop("Supplied family currently not supported. Sorry!")
       
   if(is.null(control$maxit))
      control$maxit <- 100
   if(is.null(control$tol))
      control$tol <- 1e-6
   if(is.null(control$step_halving))
      control$step_halving <- TRUE
   if(is.null(control$max_halving_steps))
      control$max_halving_steps <- 30
   if(is.null(control$trace))
      control$trace <- FALSE
   if(is.null(control$lower_limits))
      control$lower_limits <- -Inf
   if(is.null(control$upper_limits))
      control$upper_limits <- Inf

   new_coef <- cw_coef <- rep(0, num_resp*num_X)
   if(family$family %in% c("poisson","negative_binomial","tweedie")) { # Oh boy do these starting values help!!!
      sel_intercept <- seq(1, ncol(bigX),by = num_X)
      new_coef[sel_intercept] <- cw_coef[sel_intercept] <- log(colMeans(y))
      rm(sel_intercept)
      }
   if(family$family %in% c("Beta")) { 
      sel_intercept <- seq(1, ncol(bigX),by = num_X)
      new_coef[sel_intercept] <- cw_coef[sel_intercept] <- betalogitfam()$linkfun(colMeans(y))
      rm(sel_intercept)
      }
   if(!is.null(start_coef)) {
      new_coef <- cw_coef <- start_coef
      }
   new_disp_param <- cw_disp_param <- rep(0.5, num_resp)
   new_power_param <- cw_power_param <- rep(1.2, num_resp)
   
   start_fit <- list(linear_predictor = as.vector(bigX %*% cw_coef))
   start_fit$fitted <- .calc_fitted_value(start_fit$linear_predictor, family = family)

   ##-----------------------
   ## Setup working correlation matrix
   ##-----------------------
   if(corstr == "independence") 
      new_R <- cw_R <- Diagonal(n = num_resp)
   if(corstr == "unstructured") {
      get_variance <- .calc_variance(fitted = start_fit$fitted, family = family, 
                                     disp_param_long = rep(cw_disp_param, num_units), 
                                     power_param_long = rep(cw_power_param, num_units))
      new_R <- cw_R <- cov(matrix((yvec - start_fit$fitted)/sqrt(get_variance), nrow = num_units, byrow = TRUE)) # Effectively \sum_{i=1}^N (A^{-1/2}_i s_i s^T_i A^{-1/2}_i)/(N-1) 
      }
   if(corstr == "reducedrank") {
      get_variance <- .calc_variance(fitted = start_fit$fitted, family = family, 
                                    disp_param_long = rep(cw_disp_param, num_units), 
                                    power_param_long = rep(cw_power_param, num_units))
      pearsonres <- matrix((yvec - start_fit$fitted)/sqrt(get_variance), nrow = num_units, byrow = TRUE)
      rawsds <- sqrt(diag(cov(pearsonres))) 
      do_FA <- try(factanal(x = pearsonres, factors = rank, rotation = "none"), silent = TRUE)
      if(inherits(do_FA, "try-error"))
        do_FA <- try(factanal(x = pearsonres, factors = rank, rotation = "none", nstart = 100), silent = TRUE)
      new_R <- cw_R <- (tcrossprod(do_FA$loadings) + diag(x = do_FA$uniquenesses)) * (rawsds %o% rawsds)
      rm(do_FA, rawsds, pearsonres)
      }
    
   
   ##-----------------------
   ## Commence fitting
   ##-----------------------
   counter <- 1
   diff <- 1000
   alldiff <- NULL
   while(counter < control$maxit) {
      if(diff < control$tol)
         break;

      ##-----------------------
      ## Setup quantities needed to update coefficients
      ##-----------------------
      if(counter == 1) { 
         cw_W <- .calc_mu_eta(linear_predictor = start_fit$linear_predictor, family = family) ## W
         get_variance <- .calc_variance(fitted = start_fit$fitted, family = family, 
                                       disp_param_long = rep(cw_disp_param, num_units), 
                                       power_param_long = rep(cw_power_param, num_units))
         cw_Arootinv <- 1/sqrt(get_variance) ## A^{-1/2}
         rm(start_fit)
         }
      if(counter > 1) {
         cw_linear_predictor <- as.vector(bigX %*% cw_coef)
         cw_W <- .calc_mu_eta(linear_predictor = cw_linear_predictor, family = family) ## W
         cw_fitted_value <- .calc_fitted_value(cw_linear_predictor, family = family)
         get_variance <- .calc_variance(fitted = cw_fitted_value, family = family, 
                                       disp_param_long = rep(cw_disp_param, num_units), 
                                       power_param_long = rep(cw_power_param, num_units))
         cw_Arootinv <- 1/sqrt(get_variance) ## A^{-1/2}
         }
      bigRinv <- kronecker(Diagonal(n = num_units), chol2inv(chol(cw_R)))
      rm(get_variance)

      ##-----------------------
      ## Update coefficients
      ##-----------------------
      if(method == "newtonraphson") {
         get_update_coef <- .update_coef_nopen_NR(yvec = yvec, bigX = bigX, num_units = num_units, cw_coef = cw_coef, 
                                                  W = cw_W, Arootinv = cw_Arootinv, bigRinv = bigRinv, family = family)
         
         new_coef <- get_update_coef$coefficients
         temp_coef_err <- sum((new_coef - cw_coef)^2)
         if(control$step_halving & counter > 1) {
            stepsize_counter <- 1
            while(temp_coef_err > alldiff[counter-1] & stepsize_counter < control$max_halving_steps) {
               new_coef <- as.vector(cw_coef + 0.5^stepsize_counter * get_update_coef$NR_move)
               temp_coef_err <- sum((new_coef - cw_coef)^2)
               
               if(!is.finite(temp_coef_err))
                  temp_coef_err <- Inf
               
               stepsize_counter <- stepsize_counter + 1
               }
            }
         }
      if(method == "glmnet") {
         get_update_coef <- .update_coef_glmnet(cw_reparametrized_coef = cw_coef, yvec = yvec, bigXtrans = bigX, 
                                                W = cw_W, Arootinv = cw_Arootinv, Rinv = chol2inv(chol(cw_R)), 
                                                num_units = num_units, num_resp = num_resp, family = family, control = control)
         
         new_coef <- get_update_coef$coefficients
         temp_coef_err <- sum((new_coef - cw_coef)^2)
         NR_move <- new_coef - cw_coef
         if(control$step_halving & counter > 1) {
            stepsize_counter <- 1
            while(temp_coef_err > alldiff[counter-1] & stepsize_counter < control$max_halving_steps) {
               new_coef <- as.vector(cw_coef + 0.5^stepsize_counter * NR_move)
               temp_coef_err <- sum((new_coef - cw_coef)^2)
               
               if(!is.finite(temp_coef_err))
                  temp_coef_err <- Inf
               
               stepsize_counter <- stepsize_counter + 1
               }
            }
         }
   
      
      new_linear_predictor <- as.vector(bigX %*% new_coef)
      new_fitted_value <- .calc_fitted_value(linear_predictor = new_linear_predictor, family = family)
   
      
      ##-----------------------
      ## Update dispersion parameters
      ##-----------------------
      if(family$family == "gaussian") {
         rawres <- matrix(yvec - new_fitted_value, nrow = num_units, byrow = TRUE)
         new_disp_param <- as.vector(colSums(rawres^2)/(num_units - num_X)) 
         }
      if(family$family == "negative_binomial") {
         rawres <- matrix(yvec - new_fitted_value, nrow = num_units, byrow = TRUE)
         new_fitted_value_mat <- matrix(new_fitted_value, nrow = num_units, byrow = TRUE)
         for(k0 in 1:num_resp) {
            profilelogLdisp_param <- function(x) {
               varvals <- new_fitted_value_mat[,k0] + x*new_fitted_value_mat[,k0]^2
               out <- sum(-0.5*log(varvals) - 0.5*rawres[,k0]^2/varvals)
               return(out)
               }
            
            update_disp_param <- optimize(f = profilelogLdisp_param, lower = 1e-6, upper = 1e4, maximum = TRUE)
            new_disp_param[k0] <- update_disp_param$maximum
            }
         }
      if(family$family == "tweedie") {
         #new_Arootinv <- matrix(1/sqrt(new_fitted_value^rep(cw_power_param, num_units)), nrow = num_units, byrow = TRUE)
         new_fitted_value_mat <- matrix(new_fitted_value, nrow = num_units, byrow = TRUE)
         rawres <- matrix(yvec - new_fitted_value, nrow = num_units, byrow = TRUE)
         #pearsonres <- rawres * new_Arootinv
         #new_disp_param <- as.vector(colSums(pearsonres^2)/(num_units - num_X)) 

         for(k0 in 1:num_resp) {
            profilelogLpower_param <- function(x) {
               # varvals <- new_disp_param[k0]*new_fitted_value_mat[,k0]^x
               # tmp_pearsonres <- rawres[,k0]/sqrt(new_disp_param[k0]*new_fitted_value_mat[,k0]^x)
               # out <- sum(-0.5*x*log(new_fitted_value_mat[,k0]) - 0.5*tmp_pearsonres^2)
               varvals <- x[1]*new_fitted_value_mat[,k0]^x[2]
               out <- sum(-0.5*log(varvals) - 0.5*rawres[,k0]^2/varvals)
               return(out)
               }
               
            update_power_param <- optim(par = c(cw_disp_param[k0], cw_power_param[k0]), fn = profilelogLpower_param, method = "L-BFGS-B", 
                                        lower = c(1e-6, 1.001), upper = c(1e4, 2-0.001), control = list(fnscale = -1))
            new_disp_param[k0] <- update_power_param$par[1]
            new_power_param[k0] <- update_power_param$par[2]
            }
         }
      if(family$family == "Beta") {
         rawres <- matrix(yvec - new_fitted_value, nrow = num_units, byrow = TRUE)
         new_fitted_value_mat <- matrix(new_fitted_value, nrow = num_units, byrow = TRUE)
         for(k0 in 1:num_resp) {
            profilelogLdisp_param <- function(x) {
               varvals <- new_fitted_value_mat[,k0]*(1-new_fitted_value_mat[,k0])/(1+x)
               out <- sum(-0.5*log(varvals) - 0.5*rawres[,k0]^2/varvals)
               return(out)
               }
            
            update_disp_param <- optimize(f = profilelogLdisp_param, lower = 1e-6, upper = 1e4, maximum = TRUE)
            new_disp_param[k0] <- update_disp_param$maximum
            }
         }
         
      
      
      ##-----------------------
      ## Update working correlation if required
      ##-----------------------
      if(corstr == "independence") { }
      if(corstr == "unstructured") {
         get_variance <- .calc_variance(fitted = new_fitted_value, family = family, 
                                       disp_param_long = rep(new_disp_param, num_units), 
                                       power_param_long = rep(new_power_param, num_units))
         new_R <- cov(matrix((yvec - new_fitted_value)/sqrt(get_variance), nrow = num_units, byrow = TRUE)) 
         }
      if(corstr == "reducedrank") {
         get_variance <- .calc_variance(fitted = new_fitted_value, family = family, 
                                       disp_param_long = rep(new_disp_param, num_units), 
                                       power_param_long = rep(new_power_param, num_units))
         pearsonres <- matrix((yvec - new_fitted_value)/sqrt(get_variance), nrow = num_units, byrow = TRUE)
         rawsds <- sqrt(diag(cov(pearsonres))) 
         do_FA <- try(factanal(x = pearsonres, factors = rank, rotation = "none"), silent = TRUE)
         if(inherits(do_FA, "try-error"))
          do_FA <- try(factanal(x = pearsonres, factors = rank, rotation = "none", nstart = 100), silent = TRUE)
         new_R <- cw_R <- (tcrossprod(do_FA$loadings) + diag(x = do_FA$uniquenesses)) * (rawsds %o% rawsds)
         rm(do_FA, rawsds, pearsonres)
         }

        
      diff <- sum((new_coef - cw_coef)^2)
      alldiff <- c(alldiff, diff)
      if(control$trace)
         message("Iteration: ", counter, " \t Difference in estimates: ", round(diff,5))
      if(diff > 1/control$tol & counter > 10)
         break; ## Probably evidence of very poor fitting, so might as well break out and move on?
            
      cw_coef <- new_coef
      cw_disp_param <- new_disp_param
      cw_power_param <- new_power_param
      cw_R <- new_R
        
      counter <- counter + 1
      }
    
   
   ##-----------------------
   ## Fitting finished; obtaining quantities to output and tidy up. 
   ##-----------------------
   out <- list(call = match.call(), coefficients = new_coef, coefficients_matrix = matrix(new_coef, nrow = num_resp, byrow = TRUE),
               working_correlation = new_R, disp_param = cw_disp_param, power_param = cw_power_param)
   out$linear_predictor <- new_linear_predictor
   out$fitted_values <- .calc_fitted_value(linear_predictor = out$linear_predictor, family = family)
   out$residuals <- residuals_vec <- yvec - out$fitted_values
   
   cw_W <- .calc_mu_eta(linear_predictor = out$linear_predictor, family = family) ## W
   get_variance <- .calc_variance(fitted = out$fitted_values, family = family, 
                                  disp_param_long = rep(cw_disp_param, num_units), 
                                  power_param_long = rep(cw_power_param, num_units))
   cw_Arootinv <- 1/sqrt(get_variance) ## A^{-1/2}
   
   out$linear_predictor <- matrix(out$linear_predictor, nrow = num_units, byrow = TRUE)
   out$fitted_values <- matrix(out$fitted_values, nrow = num_units, byrow = TRUE)
   out$residuals <- matrix(out$residuals, nrow = num_units, byrow = TRUE)
   rownames(out$coefficients_matrix) <- colnames(out$linear_predictor) <- colnames(out$fitted_values) <- colnames(out$residuals) <- colnames(y)
   colnames(out$coefficients_matrix) <- colnames(X)
   names(out$disp_param) <- colnames(y) 
   names(out$power_param) <- colnames(y) 
   rownames(out$working_correlation) <- colnames(out$working_correlation) <- colnames(y) 
    
   
   ##-----------------------
   ## Calculating standard errors 
   ##----------------------- 
   D <- (bigX * cw_W)
   bigVinv <- Diagonal(x = cw_Arootinv) %*% kronecker(Diagonal(n = num_units), chol2inv(chol(cw_R))) %*% Diagonal(x = cw_Arootinv)
   bread_matrix_inv <- crossprod(D, bigVinv) %*% D
   bread_matrix_inv <- 0.5*(bread_matrix_inv + t(bread_matrix_inv))
   bread_matrix_inv <- chol2inv(chol(bread_matrix_inv))
   Rbar <- cov(matrix(residuals_vec * cw_Arootinv, nrow = num_units, byrow = TRUE))
   meat_matrix <- Diagonal(x = 1/cw_Arootinv) %*% kronecker(Diagonal(n = num_units), Rbar) %*% Diagonal(x = 1/cw_Arootinv)
   meat_matrix <- crossprod(D, bigVinv) %*% meat_matrix %*% (bigVinv %*% D) 
   meat_matrix <- 0.5*(meat_matrix + t(meat_matrix)) 
   sandwich_matrix <- bread_matrix_inv %*% meat_matrix %*% bread_matrix_inv 
   out$sandwich_matrix <- sandwich_matrix 
   out$stderr_matrix <- matrix(sqrt(diag(sandwich_matrix)), nrow = num_resp, byrow = TRUE) 
   rownames(out$stderr_matrix) <- colnames(y)
   colnames(out$stderr_matrix) <- colnames(X)
   
   class(out) <- "unpenalized_gee"
   return(out)
   }



   
# .manual_glm = function(bigX, y, num_resp, num_units, threshold = 1e-10, max_iter = 100) {
#    nobs <- length(y) 
#    beta = numeric(ncol(bigX))
#    beta[colSums(X) == num_units] <- log(colMeans(matrix(y, ncol = num_resp, byrow = TRUE)))
#    eta <- as.vector(bigX %*% beta)
#    mu <- poisson()$linkinv(eta)
# 
#   diff = 10000 
#   iter_count = 0
#   while(diff > threshold) {
#     W = Diagonal(x = poisson()$mu.eta(eta)) 
#     V = Diagonal(x = 1/poisson()$variance(mu))
#     good <- diag(W) != 0
#     
#     #update beta
#     beta_change = as.vector(solve(t(bigX[good,])%*%W[good,good]%*%V[good,good]%*%W[good,good]%*%bigX[good,]) %*% t(bigX[good,])%*%W[good,good]%*%V[good,good]%*%(y - mu)[good])
#     new_beta = beta + beta_change
#     # z <- eta[good] + (y - mu)[good]/poisson()$mu.eta(eta)[good]
#     # w <- (poisson()$mu.eta(eta)[good]^2)/poisson()$variance(mu)[good]
#     # fit <- lm(z ~ bigX[good,] - 1, weights = w)
#     # solve(crossprod(bigX[good,]*sqrt(w))) %*% crossprod(bigX[good,]*w, z) - fit$coefficients
#     
#     eta <- as.vector(bigX %*% new_beta)
#     mu <- poisson()$linkinv(eta)
#     beta <- new_beta
#     dev <- sum(poisson()$dev.resids(y, mu, wt = rep(1, length(y))))
#     cat("Deviance = ", dev, "\n")
# 
#     diff = sum(beta_change^2)
#     iter_count = iter_count + 1
#    }
#   
#   return(list(coefficients = beta))
#   }
# 
