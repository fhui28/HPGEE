rm(list = ls())
library(tidyverse)
library(mvtnorm)
library(mvabund)
library(glmnet)
library(cluster)  
library(ROCR)
here::i_am("simulations/setting1/setting1_binarypurepositive.R")
library(here)

source(here("code","auxfnsv0.R"))
source(here("code","calc_penaltyweights.R"))
source(here("code","create_life.R"))
source(here("code","familyfns.R"))
source(here("code","estimates.R"))
source(here("code","geefn.R"))
source(here("code","hpgeefn.R"))
source(here("code","updatecoeffns.R"))
source(here("code","ADLgeefn.R"))


num_units <- 400 ## PLEASE ALTER THIS TO THE RELEVANT NUMBER OF OBSERVATIONAL UNITS, AS DESIRED
num_resp <- 10
num_lv <- 2
num_X <- 9 # Excludes intercept
rho <- 0.5
P <- outer(1:num_X, 1:num_X, FUN = "-") %>% abs
P <- rho^P 
myfam <- binomial(link = "probit")


##---------------------
#' ## Simulate multivariate abundance data from an LVM
##---------------------
set.seed(112026)
true_lvcoefs <- matrix(runif(num_resp*num_lv, -1.5, 1.5), nrow = num_resp, ncol = num_lv)
true_disp_param <- runif(num_resp, 0, 4)
marginalization_adjustment <- matrix(sqrt(rowSums(true_lvcoefs^2) + 1), nrow = num_resp, ncol = num_X+1, byrow = FALSE)

X <- rmvnorm(num_units, sigma = P) %>% cbind(1, .)
colnames(X) <- c("Int.", paste0("x",1:num_X))


X_coefs <- cbind(runif(num_resp, -3, 0), 
                 rnorm(num_resp,0,sd=0.25), # Unique values for each response
                 replicate(3, sample(c(-0.5,0,0.5),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 3 unique values, including zero
                 replicate(3, sample(c(-1,-0.5,0.5,1),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 4 unique non-zero values
                 replicate(2, numeric(num_resp))) 
# Note no covariate is such that it has the same value for all responses -- this is not really believable in ecology =P
X_coefs <- X_coefs*marginalization_adjustment # Scaled up due to the connection between marginal and conditional regression coefficients in binary probit GEEa
X_coefs[,-1][X_coefs[,-1] < 0] <- 0
rm(P, rho, num_lv)


simy <- create_life(family = myfam, 
                    X = X, 
                    X_coefs = X_coefs, 
                    loadings = true_lvcoefs, 
                    disp_param = true_disp_param) 
true_marginal_coefs <- X_coefs/marginalization_adjustment
X <- as.data.frame(X)


##---------------------
#' ## Adaptive lasso GEEs
##---------------------
fit_adlgee <- adlgee_regpath(y = simy, 
                             X = X, 
                             family = myfam, 
                             corstr = "unstructured", 
                             multicore_path = TRUE, 
                             num_cores = detectCores() - 1)

estimates.adlgee_regpath(fit_adlgee, lambda = "BIC")$coefficients_matrix


form_lower_limits <- rep(rep(c(-Inf, 0), c(1, ncol(X)-1)), num_resp) # Basically this ensure all intercepts are unconstrained, but all covariate effects must be non-negative.
fit_adlgee_pos <- adlgee_regpath(y = simy, 
                                 X = X, 
                                 family = myfam, 
                                 corstr = "unstructured",
                                 multicore_path = TRUE, 
                                 num_cores = detectCores() - 1,
                                 control = list(lower_limits = form_lower_limits))


true_marginal_coefs
estimates.adlgee_regpath(fit_adlgee, lambda = "BIC")$coefficients_matrix
estimates.adlgee_regpath(fit_adlgee_pos, lambda = "BIC")$coefficients_matrix


data.frame(true = true_marginal_coefs[,-1] %>% as.vector,
           unconstrained_adlasso = estimates.adlgee_regpath(fit_adlgee, lambda = "BIC")$coefficients_matrix[,-1] %>% as.vector,
           constrained_lasso = estimates.adlgee_regpath(fit_adlgee_pos, lambda = "BIC")$coefficients_matrix[,-1] %>% as.vector) %>% 
   t %>% 
   dist

