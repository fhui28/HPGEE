.update_coef_nopen_NR <- function(yvec, bigX, num_units, cw_coef, W, Arootinv, bigRinv, step_size = 1, family) {
   good <- W != 0
   
   working_residual <- (yvec - .calc_fitted_value(linear_predictor = as.vector(bigX %*% cw_coef), family = family))*Arootinv
   DestArootinv <- (bigX*W*Arootinv)[good,]
   DestArootinvTbigRinv <- crossprod(DestArootinv, bigRinv[good,good]) 
   DtVinvD_inv <- solve(DestArootinvTbigRinv %*% DestArootinv)

   NR_move <- as.vector(DtVinvD_inv %*% DestArootinvTbigRinv %*% working_residual[good])
   new_coef <- as.vector(cw_coef + step_size * NR_move)

   out <- list(coefficients = as.vector(new_coef), NR_move = NR_move)
   return(out)
   }
 

.update_coef_glmnet <- function(cw_reparametrized_coef, yvec, bigXtrans, lambda = 0, penalty_weights = rep(1, length(cw_reparametrized_coef)), 
                                W, Arootinv, Rinv, num_units, num_resp, family, control) {
   bigB <- eigen(Rinv)
   bigB <- diag(x = sqrt(bigB$values)) %*% t(bigB$vectors)
   bigB <- kronecker(Diagonal(n = num_units), bigB) %*% Diagonal(x = Arootinv)

   D <- bigXtrans * W
   working_response <- D %*% cw_reparametrized_coef + (yvec - .calc_fitted_value(linear_predictor = as.vector(bigXtrans %*% cw_reparametrized_coef), family = family))
   working_response <- as.vector(working_response)
   
   #bigBD <- bigB %*% D
   #bigBz <- bigB %*% working_response
   #fit_adaptivelasso <- glmnet(x = bigBD, y = bigBz, intercept = FALSE)
   if(lambda == 0) {
      fit <- bigGlm(x = bigB %*% D, y = bigB %*% working_response, intercept = FALSE, lower.limits = control$lower_limits, upper.limits = control$upper_limits)
      new_coef <- fit$beta[,1]
      }
   if(lambda > 0) {
      fit <- glmnet(x = bigB %*% D, y = bigB %*% working_response, intercept = FALSE, lambda = lambda, penalty.factor = penalty_weights, 
                    lower.limits = control$lower_limits, upper.limits = control$upper_limits)
      new_coef <- fit$beta[,1]
      #new_coef <- coef(fit_adaptivelasso, s = lambda, exact = TRUE, x = bigB %*% D, y = bigB %*% working_response, penalty.factor = penalty_weights)[-1,1] # Using linear interpolation can speed things up per step, but often leads to longer computation time as things take more steps to complete
      }
   if(lambda == -1) { # Used purely to get glmnet to advice on a sequence of lambda
      fit <- glmnet(x = bigB %*% D, y = bigB %*% working_response, intercept = FALSE, penalty.factor = penalty_weights)
      new_coef <- fit$beta[,1]
      }
   
   out <- list(coefficients = as.vector(new_coef), fit = fit)
   return(out)
   }


# .update_coef_ncvreg <- function(cw_reparametrized_coef, yvec, bigXtrans, lambda = 0, penalty_type = "mcp", penalty_weights = rep(1, length(cw_reparametrized_coef)), 
#                                 W, Arootinv, Rinv, num_units, num_resp, family, control) {
#    bigB <- eigen(Rinv)
#    bigB <- diag(x = sqrt(bigB$values)) %*% t(bigB$vectors)
#    bigB <- kronecker(Diagonal(n = num_units), bigB) %*% Diagonal(x = Arootinv)
# 
#    D <- bigXtrans * W
#    working_response <- D %*% cw_reparametrized_coef + (yvec - .calc_fitted_value(linear_predictor = as.vector(bigXtrans %*% cw_reparametrized_coef), family = family))
#    working_response <- as.vector(working_response)
#    
#    penalty_type <- toupper(penalty_type)
#    fit_ncv <- ncvfit(X = as.matrix(bigB %*% D), y = as.vector(bigB %*% working_response), penalty = penalty_type, xtx = rep(1, ncol(bigXtrans)),  
#                      gamma = control$a, lambda = lambda, penalty.factor = penalty_weights)
#    new_coef <- fit_ncv$beta # Using linear interpolation can speed things up per step, but often leads to longer computation time as things take more steps to complete
# 
#    out <- list(coefficients = as.vector(new_coef))
#    return(out)
#    }
