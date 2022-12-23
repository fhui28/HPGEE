#' ---
#' title: Template code for running simulations involving count (Poisson) response data; see Section 4 of the associated manuscript "Simultaneous homogeneity pursuit and variable selection in regression models for multivariate abundance data" for more details 
#' abstract: Note the full simulations in the associated manuscript were run on the high performance cluster. The below template code assumes this is not available and therefore is designed for a single multi-core machine.
#' author: Originally written by FKCH
#' date: "Code started Oct 2022"
#' ---

##---------------------
#' # Load appropriate packages and files
##---------------------
rm(list = ls())
library(tidyverse)
library(mvtnorm)
library(mvabund)
library(glmnet)
library(cluster)  
library(ROCR)
here::i_am("simulations/setting1/setting1_poissonn125.R")
library(here)

source(here("code","auxfnsv0.R"))
source(here("code","calc_penaltyweights.R"))
source(here("code","create_life.R"))
source(here("code","familyfns.R"))
source(here("code","geefn.R"))
source(here("code","hpgeefn.R"))
source(here("code","updatecoeffns.R"))
source(here("code","ADLgeefn.R"))


##---------------------
#' ## Simulate multivariate abundance data from an LVM
##---------------------
num_units <- 125 ## PLEASE ALTER THIS TO THE RELEVANT NUMBER OF OBSERVATIONAL UNITS, AS DESIRED
num_resp <- 20
num_lv <- 2
num_X <- 9 # Excludes intercept
rho <- 0.5
P <- outer(1:num_X, 1:num_X, FUN = "-") %>% abs
P <- rho^P 
myfam <- poisson()

set.seed(102022)
true_lvcoefs <- matrix(runif(num_resp*num_lv, -1, 1), nrow = num_resp, ncol = num_lv)
true_disp_param <- runif(num_resp, 0, 4)

X <- rmvnorm(num_units, sigma = P) %>% cbind(1, .)
colnames(X) <- c("Int.", paste0("x",1:num_X))


X_coefs <- cbind(runif(num_resp, -2, 0), 
                 rnorm(num_resp,0,sd=0.25), # Unique values for each response
                 replicate(3, sample(c(-0.5,0,0.5),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 3 unique values, including zero
                 replicate(3, sample(c(-0.5,-0.25,0.25,0.5),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 4 unique non-zero values
                 replicate(2, numeric(num_resp))) 
# Note no covariate is such that it has the same value for all responses -- this is not really believable in ecology =P
rm(P, rho, num_lv)




dosims <- function(NAI) {
  set.seed(NAI)
          
  simy <- create_life(family = myfam, X = X, X_coefs = X_coefs, loadings = true_lvcoefs, disp_param = true_disp_param) 
  true_marginal_coefs <- X_coefs
  X <- as.data.frame(X)

     
  all_coefficient_matrices <- all_comptimes <- list()


  ##---------------------
  #' ## Proposed method: HPGEEs 
  ##---------------------
  mycorstr <- "independence"
  tic <- proc.time()
  fit_hpgee <- try(hpgee_regpath(y = simy, X = X, family = myfam, corstr = mycorstr, multicore_path = TRUE, num_cores = detectCores()-1), silent = TRUE)
  toc <- proc.time()
  if(!inherits(fit_hpgee, "try-error")) {
    fit_hpgee_BIC <- hpgee_fit(hpgee_regpath = fit_hpgee, lambda = fit_hpgee$lambda[which.min(fit_hpgee$ICs$BIC)],family = myfam, corstr = mycorstr)
    fit_hpgee_BIClogNloglogp <- hpgee_fit(hpgee_regpath = fit_hpgee, lambda = fit_hpgee$lambda[which.min(fit_hpgee$ICs$BIC_logNloglogp)], family = myfam, corstr = mycorstr)
    fit_hpgee_EBIC1 <- hpgee_fit(hpgee_regpath = fit_hpgee, lambda = fit_hpgee$lambda[which.min(fit_hpgee$ICs$EBIC1)], family = myfam, corstr = mycorstr)
    fit_hpgee_EBIC2 <- hpgee_fit(hpgee_regpath = fit_hpgee, lambda = fit_hpgee$lambda[which.min(fit_hpgee$ICs$EBIC2)], family = myfam, corstr = mycorstr)
    fit_hpgee_ERIC <- hpgee_fit(hpgee_regpath = fit_hpgee, lambda = fit_hpgee$lambda[which.min(fit_hpgee$ICs$ERIC1)], family = myfam, corstr = mycorstr)
    all_coefficient_matrices$hpgeeind_BIC <- fit_hpgee_BIC$coefficients_matrix
    all_coefficient_matrices$hpgeeind_BIClogNloglogp <- fit_hpgee_BIClogNloglogp$coefficients_matrix
    all_coefficient_matrices$hpgeeind_EBIC1 <- fit_hpgee_EBIC1$coefficients_matrix
    all_coefficient_matrices$hpgeeind_EBIC2 <- fit_hpgee_EBIC2$coefficients_matrix
    all_coefficient_matrices$hpgeeind_ERIC <- fit_hpgee_ERIC$coefficients_matrix
    
    all_comptimes$hpgeeind_BIC <- all_comptimes$hpgeeind_BIClogNloglogp <- all_comptimes$hpgeeind_EBIC1 <- all_comptimes$hpgeeind_EBIC2 <- all_comptimes$hpgeeind_ERIC <- toc-tic
    }
  if(inherits(fit_hpgee, "try-error")) {
    all_coefficient_matrices$hpgeeind_BIC <- all_coefficient_matrices$hpgeeind_BIClogNloglogp <- all_coefficient_matrices$hpgeeind_EBIC1 <- all_coefficient_matrices$hpgeeind_EBIC2 <- all_coefficient_matrices$hpgeeind_ERIC <- NA
    all_comptimes$hpgeeind_BIC <- all_comptimes$hpgeeind_BIClogNloglogp <- all_comptimes$hpgeeind_EBIC1 <- all_comptimes$hpgeeind_EBIC2 <- all_comptimes$hpgeeind_ERIC <- NA    
    }


  mycorstr <- "reducedrank"
  tic <- proc.time()
  fit_hpgee <- try(hpgee_regpath(y = simy, X = X, family = myfam, corstr = mycorstr, rank = 3, multicore_path = TRUE, num_cores = detectCores()-1), silent = TRUE)
  toc <- proc.time()
  if(!inherits(fit_hpgee, "try-error")) {
    fit_hpgee_BIC <- hpgee_fit(hpgee_regpath = fit_hpgee, lambda = fit_hpgee$lambda[which.min(fit_hpgee$ICs$BIC)], family = myfam, corstr = mycorstr, rank = 3)
    fit_hpgee_BIClogNloglogp <- hpgee_fit(hpgee_regpath = fit_hpgee, lambda = fit_hpgee$lambda[which.min(fit_hpgee$ICs$BIC_logNloglogp)], family = myfam, corstr = mycorstr, rank = 3)
    fit_hpgee_EBIC1 <- hpgee_fit(hpgee_regpath = fit_hpgee, lambda = fit_hpgee$lambda[which.min(fit_hpgee$ICs$EBIC1)], family = myfam, corstr = mycorstr, rank = 3)
    fit_hpgee_EBIC2 <- hpgee_fit(hpgee_regpath = fit_hpgee, lambda = fit_hpgee$lambda[which.min(fit_hpgee$ICs$EBIC2)], family = myfam, corstr = mycorstr, rank = 3)
    fit_hpgee_ERIC <- hpgee_fit(hpgee_regpath = fit_hpgee, lambda = fit_hpgee$lambda[which.min(fit_hpgee$ICs$ERIC1)], family = myfam, corstr = mycorstr, rank = 3)
    all_coefficient_matrices$hpgeeRR_BIC <- fit_hpgee_BIC$coefficients_matrix
    all_coefficient_matrices$hpgeeRR_BIClogNloglogp <- fit_hpgee_BIClogNloglogp$coefficients_matrix
    all_coefficient_matrices$hpgeeRR_EBIC1 <- fit_hpgee_EBIC1$coefficients_matrix
    all_coefficient_matrices$hpgeeRR_EBIC2 <- fit_hpgee_EBIC2$coefficients_matrix
    all_coefficient_matrices$hpgeeRR_ERIC <- fit_hpgee_ERIC$coefficients_matrix
    
    all_comptimes$hpgeeRR_BIC <- all_comptimes$hpgeeRR_BIClogNloglogp <- all_comptimes$hpgeeRR_EBIC1 <- all_comptimes$hpgeeRR_EBIC2 <- all_comptimes$hpgeeRR_ERIC <- toc-tic
    }
  if(inherits(fit_hpgee, "try-error")) {
    all_coefficient_matrices$hpgeeRR_BIC <- all_coefficient_matrices$hpgeeRR_BIClogNloglogp <- all_coefficient_matrices$hpgeeRR_EBIC1 <- all_coefficient_matrices$hpgeeRR_EBIC2 <- all_coefficient_matrices$hpgeeRR_ERIC <- NA
    all_comptimes$hpgeeRR_BIC <- all_comptimes$hpgeeRR_BIClogNloglogp <- all_comptimes$hpgeeRR_EBIC1 <- all_comptimes$hpgeeRR_EBIC2 <- all_comptimes$hpgeeRR_ERIC <- NA
    }
  rm(list = ls(pattern = "fit_hpgee"))
  
  
  ##---------------------
  #' ## Adaptive lasso GEEs
  ##---------------------
  mycorstr <- "independence"
  tic <- proc.time()
  fit_adlgee <- try(adlgee_regpath(y = simy, X = X, family = myfam, corstr = mycorstr, multicore_path = TRUE, num_cores = detectCores() - 1), silent = TRUE)
  toc <- proc.time()
  if(inherits(fit_adlgee, "try-error")) {
    tic <- proc.time()
    fit_adlgee <- try(adlgee_regpath(y = simy, X = X, family = myfam, corstr = mycorstr, multicore_path = FALSE, num_cores = detectCores() - 1), silent = TRUE)
    toc <- proc.time()
    }
  if(!inherits(fit_adlgee, "try-error")) {
    fit_adlgee_BIC <- adlgee_fit(adlgee_regpath = fit_adlgee, lambda = fit_adlgee$lambda[which.min(fit_adlgee$ICs$BIC)], family = myfam, corstr = mycorstr)
    fit_adlgee_BIClogNloglogp <- adlgee_fit(adlgee_regpath = fit_adlgee, lambda = fit_adlgee$lambda[which.min(fit_adlgee$ICs$BIC_logNloglogp)], family = myfam, corstr = mycorstr)
    fit_adlgee_EBIC1 <- adlgee_fit(adlgee_regpath = fit_adlgee, lambda = fit_adlgee$lambda[which.min(fit_adlgee$ICs$EBIC1)], family = myfam, corstr = mycorstr)
    fit_adlgee_EBIC2 <- adlgee_fit(adlgee_regpath = fit_adlgee, lambda = fit_adlgee$lambda[which.min(fit_adlgee$ICs$EBIC2)], family = myfam, corstr = mycorstr)
    fit_adlgee_ERIC <- adlgee_fit(adlgee_regpath = fit_adlgee, lambda = fit_adlgee$lambda[which.min(fit_adlgee$ICs$ERIC1)], family = myfam, corstr = mycorstr)
    all_coefficient_matrices$adlgeeind_BIC <- fit_adlgee_BIC$coefficients_matrix
    all_coefficient_matrices$adlgeeind_BIClogNloglogp <- fit_adlgee_BIClogNloglogp$coefficients_matrix
    all_coefficient_matrices$adlgeeind_EBIC1 <- fit_adlgee_EBIC1$coefficients_matrix
    all_coefficient_matrices$adlgeeind_EBIC2 <- fit_adlgee_EBIC2$coefficients_matrix
    all_coefficient_matrices$adlgeeind_ERIC <- fit_adlgee_ERIC$coefficients_matrix
    
    all_comptimes$adlgeeind_BIC <- all_comptimes$adlgeeind_BIClogNloglogp <- all_comptimes$adlgeeind_EBIC1 <- all_comptimes$adlgeeind_EBIC2 <- all_comptimes$adlgeeind_ERIC <- toc-tic
    }
  if(inherits(fit_adlgee, "try-error")) {
    all_coefficient_matrices$adlgeeind_BIC <- all_coefficient_matrices$adlgeeind_BIClogNloglogp <- all_coefficient_matrices$adlgeeind_EBIC1 <- all_coefficient_matrices$adlgeeind_EBIC2 <- all_coefficient_matrices$adlgeeind_ERIC <- NA
    all_comptimes$adlgeeind_BIC <- all_comptimes$adlgeeind_BIClogNloglogp <- all_comptimes$adlgeeind_EBIC1 <- all_comptimes$adlgeeind_EBIC2 <- all_comptimes$adlgeeind_ERIC <- NA
    }

  mycorstr <- "reducedrank"
  tic <- proc.time()
  fit_adlgee <- try(adlgee_regpath(y = simy, X = X, family = myfam, corstr = mycorstr, rank = 3, multicore_path = TRUE, num_cores = detectCores() - 1), silent = TRUE)
  toc <- proc.time()
  if(inherits(fit_adlgee, "try-error")) {
    tic <- proc.time()
    fit_adlgee <- try(adlgee_regpath(y = simy, X = X, family = myfam, corstr = mycorstr, rank = 3, multicore_path = FALSE, num_cores = detectCores() - 1), silent = TRUE)
    toc <- proc.time()
    }
  if(!inherits(fit_adlgee, "try-error")) {
    fit_adlgee_BIC <- adlgee_fit(adlgee_regpath = fit_adlgee, lambda = fit_adlgee$lambda[which.min(fit_adlgee$ICs$BIC)], family = myfam, corstr = mycorstr, rank = 3)
    fit_adlgee_BIClogNloglogp <- adlgee_fit(adlgee_regpath = fit_adlgee, lambda = fit_adlgee$lambda[which.min(fit_adlgee$ICs$BIC_logNloglogp)], family = myfam, corstr = mycorstr, rank = 3)
    fit_adlgee_EBIC1 <- adlgee_fit(adlgee_regpath = fit_adlgee, lambda = fit_adlgee$lambda[which.min(fit_adlgee$ICs$EBIC1)], family = myfam, corstr = mycorstr, rank = 3)
    fit_adlgee_EBIC2 <- adlgee_fit(adlgee_regpath = fit_adlgee, lambda = fit_adlgee$lambda[which.min(fit_adlgee$ICs$EBIC2)], family = myfam, corstr = mycorstr, rank = 3)
    fit_adlgee_ERIC <- adlgee_fit(adlgee_regpath = fit_adlgee, lambda = fit_adlgee$lambda[which.min(fit_adlgee$ICs$ERIC1)], family = myfam, corstr = mycorstr, rank = 3)
    all_coefficient_matrices$adlgeeRR_BIC <- fit_adlgee_BIC$coefficients_matrix
    all_coefficient_matrices$adlgeeRR_BIClogNloglogp <- fit_adlgee_BIClogNloglogp$coefficients_matrix
    all_coefficient_matrices$adlgeeRR_EBIC1 <- fit_adlgee_EBIC1$coefficients_matrix
    all_coefficient_matrices$adlgeeRR_EBIC2 <- fit_adlgee_EBIC2$coefficients_matrix
    all_coefficient_matrices$adlgeeRR_ERIC <- fit_adlgee_ERIC$coefficients_matrix

    all_comptimes$adlgeeRR_BIC <- all_comptimes$adlgeeRR_BIClogNloglogp <- all_comptimes$adlgeeRR_EBIC1 <- all_comptimes$adlgeeRR_EBIC2 <- all_comptimes$adlgeeRR_ERIC <- toc-tic
    }
  if(inherits(fit_adlgee, "try-error")) {
    all_coefficient_matrices$adlgeeRR_BIC <- all_coefficient_matrices$adlgeeRR_BIClogNloglogp <- all_coefficient_matrices$adlgeeRR_EBIC1 <- all_coefficient_matrices$adlgeeRR_EBIC2 <- all_coefficient_matrices$adlgeeRR_ERIC <- NA
    all_comptimes$adlgeeRR_BIC <- all_comptimes$adlgeeRR_BIClogNloglogp <- all_comptimes$adlgeeRR_EBIC1 <- all_comptimes$adlgeeRR_EBIC2 <- all_comptimes$adlgeeRR_ERIC <- NA
    }
  rm(list = ls(pattern = "fit_adlgee"))
  

  ##---------------------
  #' ## Separate lasso-penalized GLMs using glmnet, using 10-fold cross-validation with deviance as the loss function and selecting lambda.1se
  ##---------------------
  tic <- proc.time()
  fit_glmnet <- list(coefficient_matrix = matrix(0, nrow = num_resp, ncol = ncol(X)))
  for(k in 1:num_resp) {
    find_cv <- cv.glmnet(y = simy[,k], x = as.matrix(X[,-1,drop=FALSE]), parallel = FALSE, type.measure = "default", family = "poisson")
    fit_glmnet$coefficient_matrix[k,] <- as.vector(coef(find_cv, s = "lambda.1se", exact = TRUE))
    }
  toc <- proc.time()
  all_coefficient_matrices$glmnet <- fit_glmnet$coefficient_matrix
  all_comptimes$glmnet <- toc-tic
  rm(find_cv, fit_glmnet)


  ##---------------------
  #' ## Unpenalized GEEs
  ##---------------------
  tic <- proc.time()
  fit_geeind <- try(unpenalized_gee(y = simy, X = X, family = myfam, corstr = "independence"), silent = TRUE)
  toc <- proc.time()
  if(!inherits(fit_geeind, "try-error")) {
    all_coefficient_matrices$unpen_geeind <- fit_geeind$coefficients_matrix
    all_comptimes$unpen_geeind <- toc - tic
    }
  if(inherits(fit_geeind, "try-error")) {
    all_coefficient_matrices$unpen_geeind <- NA
    all_comptimes$unpen_geeind <- NA
    }
  tic <- proc.time()
  fit_geeRR <- try(unpenalized_gee(y = simy, X = X, family = myfam, corstr = "reducedrank", rank = 3), silent = TRUE)
  toc <- proc.time()
  if(!inherits(fit_geeRR, "try-error")) {
    all_coefficient_matrices$unpen_geeRR <- fit_geeRR$coefficients_matrix
    all_comptimes$unpen_geeRR <- toc - tic
    }
  if(inherits(fit_geeRR, "try-error")) {
    all_coefficient_matrices$unpen_geeRR <- NA
    all_comptimes$unpen_geeRR <- NA
    }

  
  ##---------------------
  #' ## Two-stage unpenalized GEEs followed by K-means clustering 
  ##---------------------
  if(!inherits(fit_geeind, "try-error")) {
    tic <- proc.time()
    fit_kmeansind <- list(coefficient_matrix = fit_geeind$coefficients_matrix)
    registerDoParallel(cores = detectCores() - 1)
    clusteringcoefficient_fn <- function(k) {
        cw_coefficients <- matrix(fit_geeind$coefficients_matrix[,k], ncol = 1) 
        gap_stat <- clusGap(cw_coefficients, FUN = kmeans, nstart = 25, K.max = num_resp-1, verbose = FALSE)
        
        cw_coefficients_cluster <- kmeans(cw_coefficients, centers = which.max(gap_stat$Tab[,3]), nstart = 25, iter.max = 50)
        out <- cw_coefficients_cluster$centers[cw_coefficients_cluster$cluster]
        return(out)
        }
    fit_kmeansind$coefficient_matrix[,-1] <- foreach(s0 = 2:ncol(X), .combine = "cbind") %dopar% clusteringcoefficient_fn(k = s0)
    toc <- proc.time()
    all_coefficient_matrices$kmeansind <- fit_kmeansind$coefficient_matrix
    all_comptimes$kmeansind <- toc-tic
    }
  if(inherits(fit_geeind, "try-error")) {
    all_coefficient_matrices$kmeansind <- NA
    all_comptimes$kmeansind <- NA
    }

  if(!inherits(fit_geeRR, "try-error")) {
    tic <- proc.time()
    fit_kmeansRR <- list(coefficient_matrix = fit_geeRR$coefficients_matrix)
    registerDoParallel(cores = detectCores() - 1)
    clusteringcoefficient_fn <- function(k) {
        cw_coefficients <- matrix(fit_geeRR$coefficients_matrix[,k], ncol = 1) 
        gap_stat <- clusGap(cw_coefficients, FUN = kmeans, nstart = 25, K.max = num_resp-1, verbose = FALSE)
        
        cw_coefficients_cluster <- kmeans(cw_coefficients, centers = which.max(gap_stat$Tab[,3]), nstart = 25, iter.max = 50)
        out <- cw_coefficients_cluster$centers[cw_coefficients_cluster$cluster]
        return(out)
        }
    fit_kmeansRR$coefficient_matrix[,-1] <- foreach(s0 = 2:ncol(X), .combine = "cbind") %dopar% clusteringcoefficient_fn(k = s0)
    toc <- proc.time()
    all_coefficient_matrices$kmeansRR <- fit_kmeansRR$coefficient_matrix
    all_comptimes$kmeansRR <- toc-tic
    }
  if(inherits(fit_geeRR, "try-error")) {
    all_coefficient_matrices$kmeansRR <- NA
    all_comptimes$kmeansRR <- NA
    }
  rm(fit_geeind, fit_geeRR, fit_kmeansind, fit_kmeansRR)
  
    
  
  ##---------------------
  #' ## Assess performance
  ##---------------------
  num_methods <- 26 # Includes oracle
  all_coefficient_matrices$true_values <- true_marginal_coefs
  true_distinct_values_nointercept <- apply(true_marginal_coefs[,-1], 2, function(x) length(unique(x))) 
  all_metrics <- matrix(NA, nrow = num_methods, ncol = 7)
  colnames(all_metrics) <- c("Sensitivity", "Specificity", "Accuracy", "F_score", "Frobenius_norm", "median_unique_values", "MAE_unique_values")
  rownames(all_metrics) <- names(all_coefficient_matrices)

  for(k0 in 1:num_methods) {
    if(is.matrix(all_coefficient_matrices[[k0]])) {
        makepred_obj <- ROCR::prediction(predictions = as.vector(1*(all_coefficient_matrices[[k0]]!=0)), labels = as.vector(1*(true_marginal_coefs != 0)))
        cw_uniquevalues_nointercept <- apply(all_coefficient_matrices[[k0]][,-1], 2, function(x) length(unique(x))) 
        
        all_metrics[k0,] <- c(
          ROCR::performance(makepred_obj, measure = "sens")@y.values[[1]][2],
          ROCR::performance(makepred_obj, measure = "spec")@y.values[[1]][2],
          ROCR::performance(makepred_obj, measure = "acc")@y.values[[1]][2], 
          ROCR::performance(makepred_obj, measure = "f")@y.values[[1]][2],
          sqrt(mean((all_coefficient_matrices[[k0]][,-1] - true_marginal_coefs[,-1])^2)),
          median(cw_uniquevalues_nointercept),
          median(abs(cw_uniquevalues_nointercept - true_distinct_values_nointercept))
          )
        }
    }
  
  
  save(all_metrics, all_comptimes, file = paste0("setting1_poissonn",num_units,"_dataset",NAI,".RData"))
  }


for(l in 1:400) 
   dosims(NAI = l)


##-------------------
sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Linux Mint 19
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
#  [1] LC_CTYPE=en_AU.UTF-8       LC_NUMERIC=C               LC_TIME=en_AU.UTF-8        LC_COLLATE=en_AU.UTF-8     LC_MONETARY=en_AU.UTF-8    LC_MESSAGES=en_AU.UTF-8   
#  [7] LC_PAPER=en_AU.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_AU.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] ROCR_1.0-7        gplots_3.0.1.2    tweedie_2.3.3     ncvreg_3.13.0     MASS_7.3-58.1     glmnet_4.1-4      Matrix_1.2-18     doParallel_1.0.17 iterators_1.0.14 
# [10] foreach_1.5.2     here_1.0.1        sf_0.9-5          blockCV_2.1.4     mgcv_1.8-40       nlme_3.1-149      mvabund_4.2.1     mvtnorm_1.1-3     forcats_0.5.2    
# [19] stringr_1.4.0     dplyr_1.0.9       purrr_0.3.4       readr_2.1.2       tidyr_1.2.0       tibble_3.1.8      ggplot2_3.3.6     tidyverse_1.3.2  
# 
# loaded via a namespace (and not attached):
#  [1] bitops_1.0-6        fs_1.5.0            lubridate_1.8.0     progress_1.2.2      httr_1.4.2          rprojroot_2.0.2     tools_3.6.3         backports_1.2.1    
#  [9] utf8_1.2.2          rgdal_1.5-16        R6_2.5.1            KernSmooth_2.23-17  DBI_1.1.0           colorspace_2.0-3    raster_3.4-5        withr_2.5.0        
# [17] sp_1.5-0            tidyselect_1.1.2    prettyunits_1.1.1   compiler_3.6.3      cli_3.4.1           rvest_1.0.3         xml2_1.3.3          labeling_0.4.2     
# [25] caTools_1.17.1.4    scales_1.2.1        classInt_0.4-3      digest_0.6.29       pkgconfig_2.0.3     dbplyr_2.2.1        rlang_1.0.6         readxl_1.3.1       
# [33] rstudioapi_0.13     shape_1.4.4         generics_0.1.3      farver_2.1.1        jsonlite_1.7.2      gtools_3.8.1        googlesheets4_1.0.1 magrittr_2.0.3     
# [41] Rcpp_1.0.8.3        munsell_0.5.0       fansi_1.0.3         abind_1.4-5         lifecycle_1.0.3     stringi_1.7.6       grid_3.6.3          gdata_2.18.0       
# [49] crayon_1.5.1        lattice_0.20-41     haven_2.5.1         splines_3.6.3       hms_1.1.2           pillar_1.8.1        codetools_0.2-16    reprex_2.0.2       
# [57] glue_1.6.2          modelr_0.1.9        vctrs_0.4.2         tzdb_0.3.0          cellranger_1.1.0    gtable_0.3.1        assertthat_0.2.1    broom_1.0.1        
# [65] e1071_1.7-3         class_7.3-17        survival_3.2-7      googledrive_2.0.0   gargle_1.2.1        units_0.6-7         statmod_1.4.36      ellipsis_0.3.2   
# 

