#' @title Calculate adaptive weights for homogeneity pursuit.
#' 
#' @param X_coefs A matrix of response-specific regression coefficients e.g., from the saturated GEE fit. The number of rows is equal to the number of responses.
#' @param kappa A vector of the length same as \code{ncol(X_coefs)} i.e., the number of covariates, which controls whether any penalization should be applied to the covariate at all. Note the default of \code{rep(c(0,1), c(1,ncol(X)-1))}, which amounts is not penalizing the first column of \code{X} i.e., the intercept column, and then penalization the other columns of \code{X} equally.
#' @param gamma1 A power parameter for the adaptive lasso penalty (for penalizing the smallest absolute coefficient). A higher value encourages more sparsity in a covariate i.e., one or more coefficients are shrunk to zero. 
#' @param gamma2 A power parameter for the adaptive lasso penalty (for penalizing successive differences in the ordered coefficient in each covariate). A higher value encourages homogeneity pursuit i.e., coefficients for a covariate are clustered into one or more groups.
#' @param gamma3 A power parameter for the adaptive weight (for penalizing the range of the coefficients in each covariate). A higher value encourages covariates that have a small range in the coefficients to be more fully homogeneous, regardless of whether that is zero-homogeneity or non-zero-homogeneity.
#' @param lambda1 A tuning parameter in the SCAD/MC+ penalty (for penalizing the smallest absolute coefficient). A higher value encourages more sparsity in a covariate i.e., one or more coefficients are shrunk to zero. 
#' @param lambda2 A tuning parameter in the SCAD/MC+ penalty (for penalizing successive differences in the ordered coefficient in each covariate). A higher value encourages homogeneity pursuit i.e., coefficients for a covariate are clustered into one or more groups.
#' @param lambda3 A tuning parameter in the SCAD/MC+ penalty (for penalizing the range of the coefficients in each covariate). A higher value encourages covariates that have a small range in the coefficients to be more fully homogeneous, regardless of whether that is zero-homogeneity or non-zero-homogeneity.
#' @param a The a parameter in the SCAD/MC+ penalty.
#' 
#' @return A matrix of adaptive weights, of the same dimensions as X_coefs. Note however that the penalty matrices apply to reparametrized coefficients.


calc_penaltyweights <- function(X_coefs, kappa = rep(1, ncol(X_coefs)),  
                                gamma1 = 1, gamma2 = 1, gamma3 = 1, 
                                lambda1 = 0, lambda2 = 0, a = 3.7) {
   
   penalty_type = "adaptive_lasso"
   #penalty_type <- match.arg(penalty_type, choices = c("adaptive_lasso", "scad", "mcp"))
# @param penalty_type The penalty type to use. Currently accepted types are (penalty functionsb based of) "adaptive_lasso" for adaptive LASSO, "scad" for SCAD, and "mcp" for MC+. 
   weights_mat <- matrix(NA, nrow = nrow(X_coefs), ncol = ncol(X_coefs))
   reparametrized_coef_mat <- X_coefs*0 
   
   for(k in 1:ncol(X_coefs)) {
      sorted_coefs <- sort(X_coefs[,k])
      reparametrised_coefs <- c(X_coefs[which.min(abs(X_coefs[,k])),k], diff(sorted_coefs))
      reparametrized_coef_mat[,k] <- reparametrised_coefs 
      
      if(penalty_type == "adaptive_lasso") {
         weights_mat[1,k] <- kappa[k]*1/abs(reparametrised_coefs[1])^gamma1
         for(k1 in 2:nrow(X_coefs))
            weights_mat[k1,k] <- kappa[k]*(1/abs(reparametrised_coefs[k1])^gamma2)*(1/abs(sum(reparametrised_coefs[2:nrow(X_coefs)]))^gamma3) # Note sum(reparametrised_coefs[2:nrow(X_coefs)]) = max(X_coefs[,k]) - min(X_coefs[,k])
         }
      
      # if(penalty_type == "scad") {
      #    if(a <= 2)
      #       stop ("The 'a' parameter in the SCAD penalty must exceed 2.")
      # 
      #    scad_derivative <- function(x, lambda) {
      #       lambda*(x <= lambda) + max(a*lambda-x, 0) / (a-1) * (x > lambda)
      #       }
      # 
      #    weights_mat[1,k] <- kappa[k]*scad_derivative(abs(reparametrised_coefs[1]), lambda = lambda1)
      #    for(k1 in 2:nrow(X_coefs))
      #       weights_mat[k1,k] <- kappa[k] * scad_derivative(abs(reparametrised_coefs[k1]) * abs(sum(reparametrised_coefs[2:nrow(X_coefs)])), lambda = lambda2) # Note sum(reparametrised_coefs[2:nrow(X_coefs)]) = max(X_coefs[,k]) - min(X_coefs[,k])
      #    }
      # 
      # if(penalty_type == "mcp") {
      #    mcp_derivative <- function(x, lambda) {
      #       lambda * max(1 - x/(a*(lambda + .Machine$double.eps)), 0)
      #       }
      # 
      #    weights_mat[1,k] <- kappa[k]*mcp_derivative(abs(reparametrised_coefs[1]), lambda = lambda1)
      #    for(k1 in 2:nrow(X_coefs))
      #       weights_mat[k1,k] <- kappa[k] * mcp_derivative(abs(reparametrised_coefs[k1]) * abs(sum(reparametrised_coefs[2:nrow(X_coefs)])), lambda = lambda2) # Note sum(reparametrised_coefs[2:nrow(X_coefs)]) = max(X_coefs[,k]) - min(X_coefs[,k])
      #    }
      
      }
   
   rownames(reparametrized_coef_mat) <- colnames(reparametrized_coef_mat) <- NULL
   return(list(penalty_weights = weights_mat, reparametrized_coefficients = reparametrized_coef_mat))
   }
