#-----------------------------
#' # Some auxilary functions
#-----------------------------
#' @title Produce a sequence of tuning parameters on the log (base 10) scale.
#' 
#' @param from Minimum tuning parameter in sequence.
#' @param to Maximum tuning parameter in sequence.
#' @param length Length of tuning parameter sequence.
#' @param decreasing Should the sequence of tuning parameters be increasing or decreasing? Defaults to \code{FALSE} i.e., increasing.
#' 
#' @return A sequence of tuning parameters.
lseq <- function (from, to, length, decreasing = FALSE) {
	stopifnot(from > 0)
	out <- 10^(seq(log10(from), log10(to), length.out = length))
	out <- out[order(out, decreasing = decreasing)]; 
	return(out)
	}


# @title Soft-threshold operator.
# 
# @param x A vector of parameters to apply the operator to.
# @param rho A vector of tuning parameters.
# 
# @return A vector of thresholded values.
.soft_threshold <- function(x, rho) {
  out <- numeric(length(x))
  
  for(k in 1:length(out)) {
    if(x[k] > 0 & rho[k] < abs(x[k])) 
          out[k] <- x[k] - rho[k] 
        if(x[k] < 0 & rho[k] < abs(x[k])) out[k] <- x[k] + rho[k] 
        if(rho[k] >= abs(x[k])) out[k] <- 0	
        }
    return(out)
    }


# @title Function that creates the reparametrization matrix
# 
# @description Basically, when the outputted matrix is right multiplied by the vector of regression coefficients (which are ordered as responses within covariates), it produces the reparametrized coefficients ready for homogeneity pursuit.
# Function adapted from and credit goes to the authors of the metafuse package
.reparametrizationmatrix <- function(X_coefs) {
   Smat <- .orderingmatrix(X_coefs = X_coefs)
   Rmat <- .adjustedifferencematrix(X_coefs = X_coefs, S = Smat)
   
   return(Rmat %*% Smat)
   }
   

# @title Function thatcreates an ordering matrix. 
# 
# @description Basically, when the outputted matrix is right multiplied by the vector of regression coefficients (which are ordered as responses within covariates), it then orders coefficients within each covariate. 
# Function adapted from and credit goes to the authors of the metafuse package
.orderingmatrix <- function(X_coefs) {
   num_resp <- nrow(X_coefs)
   num_X <- ncol(X_coefs)
   
   S_list <- list()
   for(i in 1:num_X) {
      get_order <- order(X_coefs[,i])
      Stemp <- matrix(0, nrow = num_resp, ncol = num_resp)
      for(j in 1:num_resp) {
         Stemp[j, get_order[j]] <- 1
         }
      
      S_list[[i]] <- Stemp 
      }

   out <- bdiag(S_list) 
   
   # Since out currently runs on a per-covariate basis, then reorder to run it on a per-response basis
   response_index <- rep(1:num_resp, num_X)
   out <- out[order(response_index),order(response_index)]
   return(out)
   }


# @title Function that creates an differences matrix. 
# 
# @description Basically, when the outputted matrix is right multiplied Smat, and then right multiplied by the vector of regression coefficients (which are ordered as responses within covariates), it produces the reparametrized coefficients
# Function adapted from and credit goes to the authors of the metafuse package
.adjustedifferencematrix <- function(X_coefs, S) {
   num_resp <- nrow(X_coefs)
   num_X <- ncol(X_coefs)
   new_X_coefs <- matrix(S %*% as.vector(t(X_coefs)), nrow = num_resp, byrow = TRUE)
   
   R_list <- c()
   for(i in 1:num_X) {
      Rtemp <- diag(x = c(0, rep(1, num_resp - 1)))
      Rtemp[1, which.min(abs(new_X_coefs[,i]))] <- 1
      for(k in 2:num_resp)
         Rtemp[k, k-1] <- -1
      
      R_list[[i]] <- Rtemp     
      }

   out <- bdiag(R_list)   
   
   # Since out currently runs on a per-covariate basis, then reorder to run it on a per-response basis
   response_index <- rep(1:num_resp, num_X)
   out <- out[order(response_index),order(response_index)]
   return(out)
   }
