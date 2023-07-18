#' @title Dataset generation for LVM 
#' 
#' @description This is used to generate multivariate abundance data for applications to GEE. I suppose other DGPs could be used though?
#' 
#' @param family The family/response distribution for the LVM.
#' @param X Model matrix. Note the intercept term must be manually included in X if the user wishes to include such as term.
#' @param X_coefs A matrix of response-specific regression coefficients corresponding to X. The number of rows is equal to the number of responses in the resulting simulated dataset.
#' @param true_lv An optional matrix of latent variables to use in the LVM. Otherwise, the latent variables are simulated indepndently from a standard normal distribution.
#' @param loadings A loading matrix i.e., response-specific weights, corresponing to the latent variables. The number of rows is equal to the number of responses in the resulting simulated dataset, while the number of columns is taken as the number of latent variables to use in the LVM.
#' @param offset A matrix of offsets to include, is desired. 
#' @param disp_param A response-specific vector of dispersion parameters to use for appropriate families e.g., the negative binomial family. Defaults to a vector of ones.
#' @param power_param A response-specific vector of power parameters to use for appropriate families e.g., the Tweedie family. Defaults to a vector with all values of 1.2.
#' @param max_resp The argument can be used to try and control the size of the simulated responses. Basically, for any given response the function will attempt to try and (re)simulate until all the observations are below \code{max_resp}. Currently this only appplies to Poisson, negative binomial, and Tweedie response. 

#' @return A multivariate abundance matrix of responses.


create_life <- function(family = gaussian(), X, X_coefs, lv = NULL, loadings, offset = NULL,
                        disp_param = rep(1, nrow(X_coefs)), power_param = rep(1.2, nrow(X_coefs)), max_resp = 1e8) { 
   num_units <- nrow(X)
   num_resp <- nrow(loadings)
   num_lv <- ncol(loadings)
  
   if(is.null(lv))
      lv <- matrix(rnorm(num_units*num_lv), nrow = num_units)       
   if(!is.matrix(X)) 
      X <- as.matrix(X)

   eta <- tcrossprod(X, X_coefs) + tcrossprod(lv, loadings)
   sim_y <- matrix(NA, nrow = num_units, ncol = num_resp)    
   if(!is.null(offset)) 
      eta <- eta + offset
 
 
    resp_mat <- matrix(NA, nrow = num_units, ncol = num_resp)
    for(j in 1:num_resp) {
       if(family$family == "binomial")        
          sim_y[,j] <- rbinom(num_units, size = 1, prob = pnorm(eta[,j]))
       if(family$family == "gaussian")        
          sim_y[,j] <- rnorm(num_units, mean = eta[,j], sd = sqrt(disp_param[j]))
       if(family$family == "poisson") {
          sim_y[,j] <- Inf
          try_counter <- 1
          while(any(sim_y[,j] > max_resp) & try_counter < 100) {
             sim_y[,j] <- rpois(num_units, lambda = exp(eta[,j]))
             try_counter <- try_counter + 1
             }
          }        
       if(family$family == "negative_binomial") { 
          sim_y[,j] <- Inf
          try_counter <- 1
          while(any(sim_y[,j] > max_resp) & try_counter < 100) {
             sim_y[,j] <- rnbinom(num_units, mu = exp(eta[,j]), size = 1/disp_param[j])
             try_counter <- try_counter + 1
             }
         }
       if(family$family == "tweedie") {
          sim_y[,j] <- Inf
          try_counter <- 1
          while(any(sim_y[,j] > max_resp) & try_counter < 100) {
             sim_y[,j] <- rtweedie(num_units, mu = exp(eta[,j]), phi = disp_param[j], power = power_param[j])
             try_counter <- try_counter + 1
             }
          }        
       if(family$family == "Beta")
          sim_y[,j] <- rbeta(num_units, shape1 = plogis(eta[,j])*disp_param[j], shape2 = (1-plogis(eta[,j]))*disp_param[j])
       }
    
   out <- sim_y
   colnames(out) <- paste0("resp", 1:num_resp)
   return(out)
   }

