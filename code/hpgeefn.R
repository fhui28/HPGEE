#' @title Homogeneity pursuit in GEEs
#' 
#' @description Fits a regularization path for GEEs to achieve homogeneity pursuit (clustering of coefficients within covariates) plus sparsity, assuming either an independence, an unstructured, or a reduced-rank working correlation matrix. Both homogeneity pursuit and sparsity are achieved via appropriate adapative lasso penalties, using [glmnet::glmnet()] under the hood.
#' 
#' @param y A response matrix i.e., multivariate abundance data.
#' @param X Model matrix. It is *always* assumed that the first column in \code{X} is a vector of ones corresponding to an intercept term. 
#' @param family The family to use, which then implies the link function and mean-variance relationship to use in the GEE.
#' @param nlambda The number of tuning parameters, lambda, values.
#' @param lambda_min_ratio The smallest value for the tuning parameter lambda, as a fraction of the maximum (internally-derived) lambda value (i.e., a value such that penalization leads to an intercept-only regression structure for each response). 
#' @param lambda A user-supplied tuning parameter, lambda, sequence. Note this is usually not supplied, as it is standard to have the function itself compute its own lambda sequence based on \code{nlambda} and \code{lambda_min_ratio}. If you want compute a penalized GEE for a given lambda, you can consider using \code{hpgee_fit} instead.
#' @param min_df The minimum number of non-zero coefficients in the model, including the intercept. This is useful to supply if, when the function tries to find a appropriate lambda sequence, the user wants the maximum value of lambda to not necessarily shrink to an intercept-only model. Defaults to \code{NULL}, in which case \code{min_df} is set to the number of columns of \code{y} i.e., an intercept-only model.
#' @param kappa A vector of the length same as \code{ncol(X_coefs)} i.e., the number of covariates, which controls whether any penalization should be applied to the covariate at all. Note the default of \code{rep(c(0,1), c(1,ncol(X)-1))}, which amounts is not penalizing the first column of \code{X} i.e., the intercept column, and then penalization the other columns of \code{X} equally.
#' @param corstr The structure to use for the working correlation in the GEE. Currently, the options permitted are "independence", "unstructured", and "reducedrank".
#' @param rank The rank of the reduced-rank working correlation matrix, if \code{corstr = "reducedrank"}.
#' @param num_cores The number of cores to use as  part of computation. 
#' @param multicore_path Should parallelization be used to construct the regularization path? Why this defaults to \code{FALSE}, it is often recommended to swtich this to speed things up. This is especially since construction of the regularization path of GEEs currently does *not* make use of warm starts.
#' @param control A list for controlling the GEE estimation algorithm, including the maximum number of iterations (maxit), tolerance (tol), whether to employ step halving when updating the regression coefficients (step_halving), maximum number of step_halving attempts (max_halving_steps), and whether to print updates along the way (trace).
#' 
#' @return A object of class "hpgee_regpath" with the following values, as appropriate:
#' \item{call:}{The function call.}
#' \item{lambda:}{The actual sequence of tuning parameter, lambda, values used.}
#' \item{coefficients_path:}{A sparse matrix showing the estimated regression coefficients at each point on the regularization path. The number of rows is equal \code{ncol({X}*ncol(y)}, where the coefficients run on a per-response basis, while the number of columns is equal to \code{length(lambda)}.}
#' \item{active_sets_path:}{A matrix showing the number of non-zero coefficients per covariate, at each point on the regularization path. The number of rows is equal to \code{length(lambda)}, while the number of columns is equal to \code{ncol(X)}.}
#' \item{unique_values_path:}{A matrix showing the number at each unique coefficient values per covariate, at each point on the regularization path. The number of rows is equal to \code{length(lambda)}, while the number of columns is equal to \code{ncol(X)}.}
#' \item{score_statistics_path:}{A vector showing the value of the score statistic at each point on the regularization path.}
#' \item{df_path:}{A vector showing the number of non-zero *reparameterized* coeffcients at each point on the regularization path. This is used for calculating information criterion.}
#' \item{pseudoR2:}{A vector showing the pseudo R-squared, defined as the percentage of the variation explained (as judged by the one minus the ratio of the current score statistic divided by the score statistic of the null model) at each point on the regularization path.}
#' \item{num_units:}{Equal to \code{nrow(y)}.}
#' \item{num_resp:}{Equal to \code{ncol(y)}.}
#' \item{ICs:}{A matrix of select information criteria and their values on the regularization path. The number of rows is equal to \code{length(lambda)}. The current information criteria calculated include: Akaike information criterion (AIC) with model complexity penalty of 2; Bayesian information criterion (BIC) with model complexity penalty of \code{log(num_units)};  Bayesian information criterion (BIC) with model complexity penalty of \code{log(num_units)*log(log(ncol(X)-1))}; Extended Bayesian information criterion (EBIC) with model complexity penalty of \code{log(num_units) + log(ncol(X)-1)}; Extended Bayesian information criterion (EBIC) with model complexity penalty of \code{log(num_units) + 2*log(ncol(X)-1)}; Extended regularized information criterion (ERIC) with model complexity penalty of \code{-log(lambda)}.}
#' \item{singlefit_quantities:}{A list containing quantities that can be passed into \code{hpgee_fit} to construct a single penalized GEE e.g., at a given lambda.}

function() {
   y <- simy
   nlambda = 100
   lambda_min_ratio = 1e-4
   lambda = NULL
   family = poisson()
   kappa = rep(c(0,1), c(1,ncol(X)-1))
   corstr = "independence"
   rank = 5
   num_cores = NULL
   multicore_path = FALSE
   min_df = NULL
   control = list(maxit = 100, tol = 1e-6, trace = TRUE)
   }

hpgee_regpath <- function(y, X, family = gaussian(), nlambda = 100, lambda_min_ratio = 1e-4, lambda = NULL, min_df = NULL, 
                  kappa = rep(c(0,1), c(1,ncol(X)-1)),
                  corstr = "independence", rank = 2, num_cores = NULL, multicore_path = FALSE,
                  control = list(maxit = 100, tol = 1e-6, step_halving = TRUE, max_halving_steps = 30, trace = FALSE)) {
   
   ##-----------------------
   ## Opening checks and set up
   ##-----------------------
   if(is.null(num_cores))
      num_cores <- detectCores() - 1
   registerDoParallel(cores = num_cores)
   
   num_units <- nrow(y) 
   num_resp <- ncol(y) 
   num_X <- ncol(X)
   if(is.null(colnames(y)))
      colnames(y) <- paste0("resp", 1:num_resp)
   if(is.null(colnames(X)))
      covariate_names <- paste0("x", 1:num_X)
   
   yvec <- as.vector(t(y)) # Response runs faster than unit
   X <- as.matrix(X)
   bigX <- NULL
   bigX <- foreach(i = 1:num_units, .combine = "rbind") %dopar% kronecker(Diagonal(n = num_resp),t(X[i,]))

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
   if(is.null(control$trace))
      control$trace <- FALSE
   if(is.null(control$step_halving))
      control$step_halving <- TRUE
   if(is.null(control$max_halving_steps))
      control$max_halving_steps <- 30
   if(is.null(control$lower_limits))
      control$lower_limits <- -Inf
   if(is.null(control$upper_limits))
      control$upper_limits <- Inf
   
   
   ##-----------------------
   ## Fit unpenalized GEE and calculate adaptive weightS
   ##-----------------------
   message("Fitting initial unpenalized GEE and calculating penalty weights...")
   fit_ungee <- unpenalized_gee(y = y, X = X, bigX = bigX, family = family, control = control, method = "glmnet", corstr = corstr, rank = rank, num_cores = num_cores)
   penweights_ADL <- calc_penaltyweights(X_coefs = fit_ungee$coefficients_matrix, kappa = kappa)
   penweights_ADL <- as.vector(t(penweights_ADL$penalty_weights))
   RSmatrix <- .reparametrizationmatrix(X_coefs = fit_ungee$coefficients_matrix)
   
   ##-----------------------
   ## Attempt to determine a minimum and maximum lambda, if required
   ## .hpgee_fit_manual has a built-in option to get the glmnet regularization path from the final step of the fitting algorithm. Use this as the basis to then determine a maximum lambda, and increase this until you hit an intercept-only model.
   ##-----------------------
   if(!is.null(lambda))
      lambdaseq <- lambda
   if(is.null(lambda)) {
      message("Finding an appropriate sequence of lambda values...")
      fit_togetmaxlambda <- .hpgee_fit_manual(yvec = yvec, bigX = bigX, lambda = 0, penalty_weights = penweights_ADL, 
                                          num_resp = num_resp, num_units = num_units, reparametrization_matrix = RSmatrix, getalambdaseq = TRUE,
                                          family = family, corstr = corstr, rank = rank, response_names = colnames(y), control = control)
      lambda_max <- fit_togetmaxlambda$regpath$lambda[1]
      rm(fit_togetmaxlambda)
   
      trynullfit <- .hpgee_fit_manual(yvec = yvec, bigX = bigX, lambda = lambda_max, penalty_weights = penweights_ADL, 
                           num_resp = num_resp, num_units = num_units, reparametrization_matrix = RSmatrix,
                           family = family, corstr = corstr, rank = rank, response_names = colnames(y), control = control)
      
      
      if(is.null(min_df)) 
         min_df <- num_resp
      if(min_df < num_resp)
         stop("If min_df is supplied, then it should at least be equal to the number of responses ncol(y), which corresponds to an intercept-only model.")

      while(sum(trynullfit$coefficients_matrix != 0) > min_df) {
         lambda_max <- lambda_max*1.1 
         trynullfit <- .hpgee_fit_manual(yvec = yvec, bigX = bigX, lambda = lambda_max, penalty_weights = penweights_ADL, 
                                     num_resp = num_resp, num_units = num_units, reparametrization_matrix = RSmatrix,
                                     family = family, corstr = corstr, rank = rank, response_names = colnames(y), control = control)
         }
      
      lambdaseq <- lseq(lambda_max, lambda_max*lambda_min_ratio, length = nlambda, decreasing = TRUE)
      rm(lambda_max, trynullfit)
      }
      
   ##-----------------------
   ## Construct regularization path
   ##-----------------------
   if(!multicore_path) {
      message("Constructing regularization path...")
      pb <- txtProgressBar(min = 0, max = length(lambdaseq), style = 3)
      
      coefficients_path <- Matrix(0, nrow = ncol(bigX), ncol = length(lambdaseq), sparse = TRUE)
      unique_values_path <- active_sets_path <- Matrix(0, nrow = length(lambdaseq), ncol = num_X)
      df_path <- score_statistics <- numeric(length(lambdaseq))
      AIC <- BIC <- BIC_logNloglogp <- EBIC1 <- EBIC2 <- ERIC1 <- numeric(length(lambdaseq))
      
      for(l in 1:length(lambdaseq)) {
         cwfit <- .hpgee_fit_manual(yvec = yvec, bigX = bigX, lambda = lambdaseq[l], penalty_weights = penweights_ADL, 
                                 num_resp = num_resp, num_units = num_units, reparametrization_matrix = RSmatrix,
                                 family = family, corstr = corstr, rank = rank, response_names = colnames(y), control = control)
         coefficients_path[,l,drop=FALSE] <- cwfit$coefficients 
         unique_values_path[l,,drop=FALSE] <- cwfit$unique_values
         active_sets_path[l,,drop=FALSE] <- cwfit$active_sets
         df_path[l] <- sum(cwfit$reparametrized_coefficients != 0)
         score_statistics[l] <- cwfit$score_statistic
         
         AIC[l] <- score_statistics[l] + 2*df_path[l]
         BIC[l] <- score_statistics[l] + log(num_units)*df_path[l]
         BIC_logNloglogp[l] <- score_statistics[l] + log(num_units)*log(log(num_X-1))*df_path[l] # Excludes intercept
         EBIC1[l] <- score_statistics[l] + (log(num_units) + log(num_X-1))*df_path[l]
         EBIC2[l] <- score_statistics[l] + (log(num_units) + 2*log(num_X-1))*df_path[l]
         ERIC1[l] <- score_statistics[l] - log(lambdaseq[l])*df_path[l]
         
         setTxtProgressBar(pb, l)
         }
      
      close(pb)
      rm(pb)
      }
      
   if(multicore_path) {
      message("Beginning regularization path construction (using parallel computing)...")
      
      hpgee_outerfn <- function(l) {
         cwfit <- .hpgee_fit_manual(yvec = yvec, bigX = bigX, lambda = lambdaseq[l], penalty_weights = penweights_ADL, 
                                 num_resp = num_resp, num_units = num_units, reparametrization_matrix = RSmatrix,
                                 family = family, corstr = corstr, rank = rank, response_names = colnames(y), control = control)
         
         out <- list(coefficients = cwfit$coefficients, unique_values = cwfit$unique_values, 
                     active_sets = cwfit$active_sets, df = sum(cwfit$reparametrized_coefficients != 0), score_statistic = cwfit$score_statistic)
   
         out$AIC = out$score_statistic + 2*out$df
         out$BIC = out$score_statistic + log(num_units)*out$df
         out$BIC_logNloglogp = out$score_statistic + log(num_units)*log(log(num_X-1))*out$df
         out$EBIC1 = out$score_statistic + (log(num_units) + log(num_X-1))*out$df
         out$EBIC2 = out$score_statistic + (log(num_units) + 2*log(num_X-1))*out$df
         out$ERIC1 = out$score_statistic - log(lambdaseq[l])*out$df
         
         return(out)
         }      
      
      fit_multicore_path <- foreach(l = 1:length(lambdaseq)) %dopar% hpgee_outerfn(l = l)
      
      coefficients_path <- as(sapply(fit_multicore_path, function(x) x$coefficients), "sparseMatrix")
      unique_values_path <- t(sapply(fit_multicore_path, function(x) x$unique_values))
      active_sets_path <- t(sapply(fit_multicore_path, function(x) x$active_sets))
      df_path <- sapply(fit_multicore_path, function(x) x$df)
      score_statistics <- sapply(fit_multicore_path, function(x) x$score_statistic)
      AIC <- sapply(fit_multicore_path, function(x) x$AIC) 
      BIC <- sapply(fit_multicore_path, function(x) x$BIC) 
      BIC_logNloglogp <- sapply(fit_multicore_path, function(x) x$BIC_logNloglogp) 
      EBIC1 <- sapply(fit_multicore_path, function(x) x$EBIC1) 
      EBIC2 <- sapply(fit_multicore_path, function(x) x$EBIC2) 
      ERIC1 <- sapply(fit_multicore_path, function(x) x$ERIC1)
      }
   
   
   ##-----------------------
   ## Finished fitting. Set up quantities for output 
   ##-----------------------
   out <- list(call = match.call(), lambda = lambdaseq,
               coefficients_path = coefficients_path, 
               active_sets_path = active_sets_path,
               unique_values_path = unique_values_path, 
               score_statistics_path = score_statistics, 
               df_path = df_path,
               pseudoR2 = 1 - score_statistics/score_statistics[1],
               num_units = num_units, num_resp = num_resp)
   
   out$ICs <- data.frame(AIC, BIC, BIC_logNloglogp, EBIC1, EBIC2, ERIC1)
   out$singlefit_quantities <- list(yvec = yvec, bigX = bigX, penalty_weights = penweights_ADL, 
                                    reparametrization_matrix = RSmatrix, control = control, 
                                    response_names = colnames(y), covariate_names = colnames(X))
   
   class(out) <- "hpgee_regpath"
   return(out)
   }



#' @title Homogeneity pursuit in GEEs (single fit)
#' 
#' @description Fits single penalized GEEs to achieve homogeneity pursuit (clustering of coefficients within covariates) plus sparsity, assuming either an independence, an unstructured, or a reduced-rank working correlation matrix. Both homogeneity pursuit and sparsity are achieved via appropriate adapative lasso penalties, using [glmnet::glmnet()] under the hood. 
#' This function is designed to used after the regularization path has been constructed using the \code{hpgee_regpath} function. That is, after applying \code{hpgee_regpath}, this function is then used to calculate the GEE fit at specific values of the tuning parameter, lambda.
#' 
#' @param hpgee_regpath An object of class "hpgee_path".
#' @param lambda A user-supplied tuning parameter, lambda. Only a scalar value is permitted here.
#' @param family The family to use, which then implies the link function and mean-variance relationship to use in the GEE.
#' @param corstr The structure to use for the working correlation in the GEE. Currently, the options permitted are "independence", "unstructured", and "reducedrank".
#' @param rank The rank of the reduced-rank working correlation matrix, if \code{corstr = "reducedrank"}.
#' @param num_cores The number of cores to use as part of computation. 
#' 
#' @return A object of class "hpgee_fit" with the following values, as appropriate:
#' \item{call:}{The function call.}
#' \item{coefficients:}{The full vector of regression coefficients (running on a per-response basis).}
#' \item{coefficients_matrix:}{The regression coefficients but arranged into a matrix, with the number of rows equal to the number of responses i.e., \code{ncol(y)}.}
#' \item{reparametrized_coefficients:}{The full vector of estimated *reparametrized* regression coefficients. This can be safely ignored for practical purposes. }
#' \item{working_correlation:}{The working correlation matrix.}
#' \item{disp_param:}{The response-specific dispersion parameters. Note families which do not require dispersion parameters will still return a vector of values, although these are not updated/can be ignored.}
#' \item{power_param:}{The response-specific power parameters. Note families which do not require power parameters will still return a vector of values, although these are not updated/can be ignored.}
#' \item{active_sets:}{A vector showing the number of non-zero coefficients per covariate.}
#' \item{unique_values:}{A vector showing the number of unique coefficients per covariate.}
#' \item{linear_predictor:}{The matrix of linear predictors. The number of columns is equal to the number of responses i.e., \code{ncol(y)}.}
#' \item{fitted_values:}{The matrix of fitted values. The number of columns is equal to the number of responses i.e., \code{ncol(y)}.}
#' \item{residuals:}{The matrix of raw residuals. The number of columns is equal to the number of responses i.e., \code{ncol(y)}.}


## TODO: Standard errors for beta. This is quite difficult though: oracle allows you to get standard errors for the non-zero reparametrized coefficients. But you need the joint distribution of both zero and non-zero reparametrized coefficients to get the joint distribution of the original coefficients. Bootstrap is not expected to work well, and we already know from the Poestcher and Leeb work that finite sample distributions are pretty horrible anyway!

hpgee_fit <- function(hpgee_regpath, lambda, family = gaussian(), corstr = "independence", rank = 2) {
   if(length(lambda) > 1)
      stop("Only a scalar value of the tuning parameter, lambda, is permitted.")
   
   out <- .hpgee_fit_manual(yvec = hpgee_regpath$singlefit_quantities$yvec, 
                            bigX = hpgee_regpath$singlefit_quantities$bigX,
                            lambda = lambda,
                            penalty_weights = hpgee_regpath$singlefit_quantities$penalty_weights,
                            num_resp = hpgee_regpath$num_resp,
                            num_units = hpgee_regpath$num_units,
                            reparametrization_matrix = hpgee_regpath$singlefit_quantities$reparametrization_matrix,
                            family = family, 
                            corstr = corstr,
                            rank = rank,
                            response_names = hpgee_regpath$singlefit_quantitles$response_names,
                            covariate_names = hpgee_regpath$singlefit_quantitles$covariate_names,
                            control = hpgee_regpath$singlefit_quantities$control)
   
   out$score_statistic <- NULL
   return(out)
   }


function() {
   lambda = 0
   penalty_weights = penweights_ADL
   reparametrization_matrix = RSmatrix
   family = gaussian()
   lambda = 0 
   control = list(trace = TRUE)
   start_coef = NULL
   response_names = colnames(y)
   covariate_names <- NULL
   getalambdaseq <- TRUE
   }

#' For yvec, it is assumed that response runs faster than site
.hpgee_fit_manual <- function(yvec, bigX, lambda = 0, penalty_weights, 
                          num_resp, num_units, reparametrization_matrix, 
                          family = gaussian(), corstr = "independence", rank = 2, start_coef = NULL, 
                          response_names = NULL, covariate_names = NULL, getalambdaseq = FALSE, control) {
   
   ##-----------------------
   ## Opening checks and set up
   ##-----------------------
   bigXtrans <- bigX %*% solve(reparametrization_matrix)
   num_X <- ncol(bigXtrans)/num_resp

   penalty_type <- "adaptive_lasso"
   #penalty_type <- match.arg(penalty_type, choices = c("adaptive_lasso","scad", "mcp")) 

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
   #if(is.null(control$a))
   #   control$a <- 2
   if(is.null(control$trace))
      control$trace <- FALSE
   if(is.null(control$step_halving))
      control$step_halving <- TRUE
   if(is.null(control$max_halving_steps))
      control$max_halving_steps <- 30
   if(is.null(control$lower_limits))
      control$lower_limits <- -Inf
   if(is.null(control$upper_limits))
      control$upper_limits <- Inf
   
   if(is.null(response_names))
      response_names <- paste0("resp", 1:num_resp)
   if(is.null(covariate_names))
      covariate_names <- paste0("x", 1:num_X)
   
   new_coef_mat <- cw_coef_mat <- matrix(0, nrow = num_resp, ncol = num_X)
   if(family$family %in% c("poisson","negative_binomial","tweedie")) { # Oh boy do these starting values help!!!
      new_coef_mat[,1] <- cw_coef_mat[,1] <- log(colMeans(matrix(yvec, nrow = num_units, byrow = TRUE)))
      }
   if(family$family %in% c("Beta")) { 
      new_coef_mat[,1] <- cw_coef_mat[,1] <- betalogitfam()$linkfun(colMeans(matrix(yvec, nrow = num_units, byrow = TRUE)))
      }
   if(!is.null(start_coef)) {
      new_coef_mat <- cw_coef_mat <- matrix(start_coef, nrow = num_resp)
      }
   new_reparametrized_coef <- cw_reparametrized_coef <- (reparametrization_matrix %*% as.vector(t(new_coef_mat)))
   new_disp_param <- cw_disp_param <- rep(0.5, num_resp)
   new_power_param <- cw_power_param <- rep(1.2, num_resp)
   
   start_fit <- list(linear_predictor = as.vector(bigXtrans %*% cw_reparametrized_coef))
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
         cw_linear_predictor <- as.vector(bigXtrans %*% cw_reparametrized_coef)
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
      if(penalty_type == "adaptive_lasso") {
         get_update_coef <- .update_coef_glmnet(cw_reparametrized_coef = cw_reparametrized_coef, yvec = yvec, bigXtrans = bigXtrans, 
                                                lambda = lambda, penalty_weights = penalty_weights,
                                                W = cw_W, Arootinv = cw_Arootinv, Rinv = chol2inv(chol(cw_R)), 
                                                num_units = num_units, num_resp = num_resp, family = family, control = control)
         
         new_reparametrized_coef <- get_update_coef$coefficients
         temp_coef_err <- sum((new_reparametrized_coef - cw_reparametrized_coef)^2)
         NR_move <- new_reparametrized_coef  - cw_reparametrized_coef
         if(control$step_halving & counter > 1) {
            stepsize_counter <- 1
            while(temp_coef_err > alldiff[counter-1] & stepsize_counter < control$max_halving_steps) {
               new_reparametrized_coef <- as.vector(cw_reparametrized_coef + 0.5^stepsize_counter * NR_move)
               temp_coef_err <- sum((new_reparametrized_coef - cw_reparametrized_coef)^2)
               
               if(!is.finite(temp_coef_err))
                  temp_coef_err <- Inf
               
               stepsize_counter <- stepsize_counter + 1
               }
            }
         }
      # if(penalty_type %in% c("mcp","scad")) {
      #    get_update_coef <- .update_coef_ncvreg(cw_reparametrized_coef = cw_reparametrized_coef, yvec = yvec, bigXtrans = bigXtrans, 
      #                                           lambda = lambda, penalty_type = penalty_type, penalty_weights = penalty_weights,
      #                                           W = cw_W, Arootinv = cw_Arootinv, Rinv = chol2inv(chol(cw_R)), 
      #                                           num_units = num_units, num_resp = num_resp, family = family, control = control)
      #    }
      new_coef_mat <- matrix(solve(reparametrization_matrix) %*% new_reparametrized_coef, nrow = num_resp, byrow = TRUE)

      new_linear_predictor <- as.vector(bigXtrans %*% new_reparametrized_coef)
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
         if(inherits(do_FA, "try-error"))
          do_FA <- try(factanal(x = pearsonres, factors = rank, rotation = "none", nstart = 100, lower = 0.01), silent = TRUE)
         if(inherits(do_FA, "try-error"))
          do_FA <- try(factanal(x = pearsonres, factors = rank, rotation = "none", nstart = 100, lower = 0.1), silent = TRUE)
         new_R <- cw_R <- (tcrossprod(do_FA$loadings) + diag(x = do_FA$uniquenesses)) * (rawsds %o% rawsds)
         rm(do_FA, rawsds, pearsonres)
         }

        
      ##-----------------------
      ## Finish up iteration
      ##-----------------------
      diff <- sum((new_coef_mat - cw_coef_mat)^2) # Checking convergence on the original coefficients
      alldiff <- c(alldiff, diff)
      if(control$trace)
         message("Iteration: ", counter, " \t Difference in (reparametrized) estimates: ", round(diff,4))
      if(diff > 1/control$tol & counter > 10)
         break; ## Probably evidence of very poor fitting, so might as well break out and move on?
            
      cw_reparametrized_coef <- new_reparametrized_coef
      cw_coef_mat <- new_coef_mat
      cw_disp_param <- new_disp_param
      cw_power_param <- new_power_param
      cw_R <- new_R
        
      counter <- counter + 1
      }
    
   
   ##-----------------------
   ## Fitting finished; obtaining quantities to output and tidy up. 
   ##-----------------------
   out <- list(call = match.call(), coefficients = as.vector(t(cw_coef_mat)), coefficients_matrix = cw_coef_mat, 
               reparametrized_coefficients = cw_reparametrized_coef,
               working_correlation = cw_R, disp_param = cw_disp_param, power_param = cw_power_param)
   
   out$active_sets <- apply(new_coef_mat != 0, 2, sum)
   out$unique_values <- apply(new_coef_mat, 2, function(x) length(unique(x)))
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
   
   names(out$active_sets) <- names(out$unique_values) <- covariate_names
   rownames(out$coefficients_matrix) <- colnames(out$linear_predictor) <- colnames(out$fitted_values) <- colnames(out$residuals) <- response_names
   colnames(out$coefficients_matrix) <- covariate_names
   names(out$disp_param) <- response_names 
   names(out$power_param) <- response_names 
   rownames(out$working_correlation) <- colnames(out$working_correlation) <- response_names 
   
   if(getalambdaseq) {
      formregpath <- .update_coef_glmnet(cw_reparametrized_coef = cw_reparametrized_coef, yvec = yvec, bigXtrans = bigXtrans, 
                                         lambda = -1, penalty_weights = penalty_weights,
                                         W = cw_W, Arootinv = cw_Arootinv, Rinv = chol2inv(chol(cw_R)), 
                                         num_units = num_units, num_resp = num_resp, family = family, control = control)
      out$regpath <- formregpath$fit   
      }

   
   ##-----------------------
   ## Calculate scores and thus a score statistic
   ##-----------------------
   D <- bigXtrans * cw_W
   bigVinv <- Diagonal(x = cw_Arootinv) %*% kronecker(Diagonal(n = num_units), chol2inv(chol(cw_R))) %*% Diagonal(x = cw_Arootinv)
   score_vec <- crossprod(D, bigVinv %*% residuals_vec)
   Rbar <- cov(matrix(residuals_vec * cw_Arootinv, nrow = num_units, byrow = TRUE))
   score_covariance <- Diagonal(x = 1/cw_Arootinv) %*% kronecker(Diagonal(n = num_units), Rbar) %*% Diagonal(x = 1/cw_Arootinv)
   score_covariance <- crossprod(D, bigVinv) %*% score_covariance %*% (bigVinv %*% D) 
   score_covariance <- 0.5*(score_covariance + t(score_covariance)) 
   score_covariance_inv <- try(chol2inv(chol(score_covariance)), silent = TRUE)
   if(inherits(score_covariance_inv, "try-error")) {
       score_covariance <- Matrix::nearPD(score_covariance)$mat
       score_covariance_inv <- chol2inv(chol(score_covariance))
       }
   out$score_statistic <- as.vector(crossprod(score_vec, score_covariance_inv) %*% score_vec)
   rm(D, bigVinv, score_vec, Rbar, score_covariance)
   
   
   class(out) <- "hpgee_fit"
   return(out)
   }

 
