.calc_fitted_value <- function(linear_predictor, family) {
   if(family$family == "gaussian")
      out <- linear_predictor
   if(family$family %in% c("poisson","negative_binomial","tweedie"))
      out <- poisson()$linkinv(linear_predictor)
   if(family$family == "binomial")
      out <- family$linkinv(linear_predictor)
   if(family$family == "Beta")
      out <- family$linkinv(linear_predictor)
   
   return(out)
   }


.calc_variance <- function(fitted, family, disp_param_long, power_param_long) {
   if(family$family == "gaussian")
      out <- disp_param_long
   if(family$family %in% c("poisson"))
      out <- family$variance(fitted)
   if(family$family %in% c("negative_binomial"))
      out <- family$variance(mu = fitted, phi = disp_param_long)
   if(family$family %in% c("tweedie"))
      out <- family$variance(mu = fitted, power = power_param_long, phi = disp_param_long)
   if(family$family == "binomial")
      out <- family$variance(fitted)
   if(family$family == "Beta")
      out <- family$variance(mu = fitted, phi = disp_param_long)
   
   return(out)
   }


.calc_mu_eta <- function(linear_predictor, family) {
   if(family$family == "gaussian")
      out <- family$mu.eta(linear_predictor)
   if(family$family %in% c("poisson","negative_binomial","tweedie"))
      out <- poisson()$mu.eta(linear_predictor)
   if(family$family == "binomial")
      out <- family$mu.eta(linear_predictor)
   if(family$family == "Beta")
      out <- binomial(link = "logit")$mu.eta(linear_predictor)
   
   return(out)
   }




#' @title Family object for beta distribution.
#' 
#' @description A specified family object for the beta distribution using the logit link function. 
#'
#' @return An object of class "family".
betalogitfam <- function() {
   link <- "logit"
   linkfun <- function(mu) 
      return(log(mu/(1-mu)))
   linkinv <- binomial()$linkinv 
   mu.eta <- binomial()$mu.eta 
   variance <- function(mu, phi) 
      return(mu*(1-mu)/(1+phi))

   structure(list(family = "Beta", link = link, linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance), class = "family")
   }


#' @title Family object for negative binomial distribution.
#' 
#' @description  A specified family object for the negative binomial distribution using the log link function and a quadratic mean-variance relationship.
#'
#' @return An object of class "family".
nb2 <- function() {
	link <- "log"
	linkfun <- function(mu) 
	   return(log(mu))
	linkinv <- function(eta) 
	   return(pmax(exp(eta), .Machine$double.eps))
	mu.eta <- function(eta) 
	   return(pmax(exp(eta), .Machine$double.eps))
	variance <- function(mu, phi) 
	   return(pmax(mu+phi*mu^2, .Machine$double.eps))
  
	structure(list(family = "negative_binomial", link = link, linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance), class = "family")
	}
	
	
#' @title Family object for Tweedie distribution.
#' 
#' @description A specified family object for the Tweedie distribution using the log link function and a power mean-variance relationship.
#'
#' @return An object of class "family". 
tweedielogfam <- function() {
   link <- "log"
   linkfun <- function(mu) 
      return(log(mu))
   linkinv <- function(eta) 
      return(pmax(exp(eta), .Machine$double.eps))
   mu.eta <- function(eta) 
      return(pmax(exp(eta), .Machine$double.eps))
   variance <- function(mu, power, phi) 
      return(pmax(phi * mu^power, .Machine$double.eps))

   structure(list(family = "tweedie", link = link, linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance), class = "family")
   }

