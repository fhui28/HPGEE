#' @title Extract estimates from a homogeneity pursuit GEE regularization path.
#'
#' @description
#' Extracts the estimated regression coefficients, along with the number of non-zero coefficients (active sets) and the number of unique coefficient values per covariate, from an object of class `hpgee_regpath` at a chosen value of the tuning parameter.
#'
#' @param object An object of class `hpgee_regpath`, as produced by [hpgee_regpath()].
#' @param lambda The tuning parameter at which estimates are extracted. This can either be a string corresponding to one of the information criteria stored in `object$ICs`, in which case the point on the regularization path minimizing that criterion is used, or a single numeric value. Similar to [glmnet::coef.glmnet()], in the case of the latter linear interpolation is used to construct estimates for values of `lambda` that do not coincide with those used in constructing the regularization path. Defaults to \code{"BIC"}.
#' @param ... Not used.
#'
#' @return A list with the following elements:
#' \item{coefficients_matrix}{The matrix of estimated regression coefficients, with the number of rows equal to the number of responses and the number of columns equal to the number of covariates (including the intercept).}
#' \item{coefficients}{The full vector of estimated regression coefficients, running on a per-response basis (consistent with the rows of \code{object$coefficients_path}).}
#' \item{active_sets}{A vector showing the number of non-zero coefficients per covariate.}
#' \item{unique_values}{A vector showing the number of unique coefficient values per covariate.}
#' \item{lambda}{The numerical value of the tuning parameter corresponding to the extracted estimates.}
#'

estimates.hpgee_regpath <- function(object, lambda = "BIC", ...) {
   if(!inherits(object, "hpgee_regpath"))
      stop("object must be of class 'hpgee_regpath'.")
   
   num_resp <- object$num_resp
   num_X <- ncol(object$active_sets_path)
   
   response_names <- object$singlefit_quantities$response_names
   if(is.null(response_names))
      response_names <- paste0("resp", 1:num_resp)
   covariate_names <- object$singlefit_quantities$covariate_names
   if(is.null(covariate_names))
      covariate_names <- paste0("x", 1:num_X)
   
   
   ##-------------------
   ## Extract coefficients
   ##-------------------
   if(!is.numeric(lambda)) {
      lambda <- match.arg(lambda, choices = colnames(object$ICs))
      sel_lambda <- which.min(object$ICs[, lambda])
      
      coefficients <- as.vector(object$coefficients_path[, sel_lambda])
      active_sets <- as.vector(object$active_sets_path[sel_lambda, ])
      unique_values <- as.vector(object$unique_values_path[sel_lambda, ])
      final_lambda <- object$lambda[sel_lambda]
      }
   
   if(is.numeric(lambda)) {
      if(length(lambda) != 1)
         stop("Only a scalar value of the tuning parameter, lambda, is permitted.")
      
      lambda_list <- .lambda_interp(lambda_seq = object$lambda, s = lambda)
      coefficients <- as.vector(object$coefficients_path[, lambda_list$left] * lambda_list$frac + object$coefficients_path[, lambda_list$right] * (1 - lambda_list$frac))
      final_lambda <- min(max(lambda, min(object$lambda)), max(object$lambda)) # Values outside the path are clamped to its ends
      }
   
   
   ##-------------------
   ## Format
   ##-------------------
   coefficients_matrix <- matrix(coefficients, nrow = num_resp, byrow = TRUE)
   rownames(coefficients_matrix) <- response_names
   colnames(coefficients_matrix) <- covariate_names
   names(coefficients) <- paste0(rep(response_names, each = num_X), ":", rep(covariate_names, num_resp))
   
   if(is.numeric(lambda)) {
      active_sets <- colSums(coefficients_matrix != 0)
      unique_values <- apply(coefficients_matrix, 2, function(x) length(unique(x)))
      }
   names(active_sets) <- names(unique_values) <- covariate_names
   
   out <- list(coefficients_matrix = coefficients_matrix,
               coefficients = coefficients,
               active_sets = active_sets,
               unique_values = unique_values,
               lambda = final_lambda)
   return(out)
   }


## Linear interpolation along a (decreasing) lambda sequence; adapted from the internal lambda.interp function in glmnet.
## Estimates at s are then constructed as path[, left] * frac + path[, right] * (1 - frac).
.lambda_interp <- function(lambda_seq, s) {
   if(length(lambda_seq) == 1){
      nums <- length(s)
      left <- rep(1,nums)
      right <- left
      sfrac <-rep(1,nums)
      } else {
      k <- length(lambda_seq)
      sfrac <- (lambda_seq[1] - s)/(lambda_seq[1] - lambda_seq[k])
      lambda_seq <- (lambda_seq[1] - lambda_seq)/(lambda_seq[1] - lambda_seq[k])
      
      sfrac[sfrac < min(lambda_seq)] <- min(lambda_seq)
      sfrac[sfrac > max(lambda_seq)] <- max(lambda_seq)
      
      coord <- stats::approx(lambda_seq, seq(lambda_seq), sfrac)$y
      left <- floor(coord)
      right <- ceiling(coord)
      
      sfrac <- (sfrac - lambda_seq[right])/(lambda_seq[left] - lambda_seq[right])
      sfrac[left == right] <- 1
      sfrac[abs(lambda_seq[left] - lambda_seq[right]) < .Machine$double.eps] <- 1
      }
   
   return(list(left = left,right = right,frac = sfrac))
   }
