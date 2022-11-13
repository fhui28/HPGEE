#' ---
#' title: Code to process results for simulation study; see Section 4 of the associated manuscript "Simultaneous homogeneity pursuit and variable selection in regression models for multivariate abundance data" for more details
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
here::i_am("simulations/setting1/processresults.R")
library(here)


##---------------------
#' ## Gather results
##---------------------
num_resp <- 20
num_lv <- 2
num_X <- 9 # Excludes intercept
rho <- 0.5
P <- outer(1:num_X, 1:num_X, FUN = "-") %>% abs
P <- rho^P 
myfam <- "poisson" ## PLEASE CHANGE THIS TO THE RELEVANT RESPONSE TYPE


num_units <- c(125,250,500,1000)
num_dats <- 400
num_methods <- 26 # Includes oracle
combined_results <- array(NA, dim = c(length(num_units), num_dats, num_methods, 7))
combined_comptime <- array(NA, dim = c(length(num_units), num_dats, num_methods-1))

for(l0 in 1:length(num_units)) { for(k0 in 1:num_dats) {
   filename <- paste0(myfam, "setting1_", myfam, "n", num_units[l0], "_dataset", k0, ".RData")
   
   if(file.exists(filename)) {
      message("Gathering results from file ", filename)
      load(file = filename)      
      
      combined_results[l0,k0,,] <- all_metrics
      combined_comptime[l0,k0,] <- sapply(all_comptimes, function(x) x[3])
      
      if(l0 == 1 & k0 == 1) {
         dimnames(combined_results) <- list(N = num_units, dataset = 1:num_dats, methods = rownames(all_metrics), metrics = colnames(all_metrics))
         dimnames(combined_comptime) <- list(N = num_units, dataset = 1:num_dats, names(all_comptimes))
         }
      
      rm(all_metrics, all_comptimes)
      }
   } }



##---------------------
#' ## Explore and plot results
##---------------------
apply(combined_results, c(4,3,1), mean, na.rm = TRUE) %>% 
   ftable %>% 
   round(3)


apply(combined_results, c(4,3,1), sd, na.rm = TRUE) %>% 
   ftable %>% 
   round(3)



combined_results_long <- combined_results %>%
   as.data.frame.table
colnames(combined_results_long) <- c("N", "Dataset", "Method", "Metric", "Value")
combined_results_long <- combined_results_long %>% 
   mutate(N = fct_inorder(N), Method = fct_inorder(Method), Metric = fct_inorder(Metric))  
combined_results_long$Method <- fct_recode(combined_results_long$Method,
                                            "HPGEE (Ind. and log(N))" = "hpgeeind_BIC", 
                                            "HPGEE (Ind. and log(N)loglog(p))" = "hpgeeind_BIClogNloglogp", 
                                            "HPGEE (Ind. and log(N) + log(p))" = "hpgeeind_EBIC1", 
                                            "HPGEE (Ind. and log(N) + 2log(p))" = "hpgeeind_EBIC2", 
                                            "HPGEE (Ind. and -log(lambda))" = "hpgeeind_ERIC", 
                                            
                                            "HPGEE (RR and log(N))" = "hpgeeRR_BIC", 
                                            "HPGEE (RR and log(N)loglog(p))" = "hpgeeRR_BIClogNloglogp", 
                                            "HPGEE (RR and log(N) + log(p))" = "hpgeeRR_EBIC1", 
                                            "HPGEE (RR and log(N) + 2log(p))" = "hpgeeRR_EBIC2", 
                                            "HPGEE (RR and -log(lambda))" = "hpgeeRR_ERIC", 
                                            
                                            "ADLGEE (Ind. and log(N))" = "adlgeeind_BIC", 
                                            "ADLGEE (Ind. and log(N)loglog(p))" = "adlgeeind_BIClogNloglogp", 
                                            "ADLGEE (Ind. and log(N) + log(p))" = "adlgeeind_EBIC1", 
                                            "ADLGEE (Ind. and log(N) + 2log(p))" = "adlgeeind_EBIC2", 
                                            "ADLGEE (Ind. and -log(lambda))" = "adlgeeind_ERIC", 
                                            
                                            "ADLGEE (RR and log(N))" = "adlgeeRR_BIC", 
                                            "ADLGEE (RR and log(N)loglog(p))" = "adlgeeRR_BIClogNloglogp", 
                                            "ADLGEE (RR and log(N) + log(p))" = "adlgeeRR_EBIC1", 
                                            "ADLGEE (RR and log(N) + 2log(p))" = "adlgeeRR_EBIC2", 
                                            "ADLGEE (RR and -log(lambda))" = "adlgeeRR_ERIC", 
                                            
                                            "glmnet" = "glmnet",
                                            
                                            "Unpen. GEE (Ind.)" = "unpen_geeind", 
                                            "Unpen. GEE (RR)" = "unpen_geeRR", 
                                            
                                            "K-means (Ind.)" = "kmeansind", 
                                            "K-means (RR)" = "kmeansRR") 
combined_results_long$N <- fct_recode(combined_results_long$N,
                                      "N = 125" = "125",
                                      "N = 250" = "250",
                                      "N = 500" = "500",
                                      "N = 1000" = "1000")
combined_results_long$methods_group <- combined_results_long$Method %>% 
   fct_collapse(HPGEE = levels(combined_results_long$Method)[1:10],
                ADLGEE = levels(combined_results_long$Method)[11:20],
                glmnet = "glmnet",
                UPGEE = levels(combined_results_long$Method)[22:23],
                Kmeans = levels(combined_results_long$Method)[24:25])



ggplot(combined_results_long %>% 
          dplyr::filter(Metric == "Sensitivity" & Method != "true_values") %>% 
          group_by(N, Method, Metric, methods_group) %>% 
          summarise(mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE)), 
       aes(x = Method, y = mean, ymin = mean-sd, ymax = mean+sd, color = methods_group)) +
   geom_pointrange(show.legend = FALSE) +
   scale_color_viridis_d() +
   facet_wrap(N ~ ., ncol = 2, scales = "free_x") +
   coord_flip() +
   theme_bw() +
   labs(y = "Sensitivity")


ggplot(combined_results_long %>% 
          dplyr::filter(Metric == "Specificity" & Method != "true_values") %>% 
          group_by(N, Method, Metric, methods_group) %>% 
          summarise(mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE)), 
       aes(x = Method, y = mean, ymin = mean-sd, ymax = mean+sd, color = methods_group)) +
   geom_pointrange(show.legend = FALSE) +
   scale_color_viridis_d() +
   facet_wrap(N ~ ., ncol = 2, scales = "free_x") +
   coord_flip() +
   theme_bw() +
   labs(y = "Specificity")


ggplot(combined_results_long %>% 
          dplyr::filter(Metric == "Accuracy" & Method != "true_values") %>%  
          group_by(N, Method, Metric, methods_group) %>% 
          summarise(mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE)), 
       aes(x = Method, y = mean, ymin = mean-sd, ymax = mean+sd, color = methods_group)) +
   geom_pointrange(show.legend = FALSE) +
   scale_color_viridis_d() +
   facet_wrap(N ~ ., ncol = 2, scales = "free_x") +
   coord_flip() +
   theme_bw() +
   labs(y = "Accuracy")


ggplot(combined_results_long %>% 
          dplyr::filter(Metric == "F_score" & Method != "true_values") %>% 
          group_by(N, Method, Metric, methods_group) %>% 
          summarise(mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE)), 
       aes(x = Method, y = mean, ymin = mean-sd, ymax = mean+sd, color = methods_group)) +
   geom_pointrange(show.legend = FALSE) +
   scale_color_viridis_d() +
   facet_wrap(N ~ ., ncol = 2, scales = "free_x") +
   coord_flip() +
   theme_bw() +
   labs(y = "F-score")


ggplot(combined_results_long %>% 
          dplyr::filter(Metric == "Frobenius_norm" & Method != "true_values") %>%  
          group_by(N, Method, Metric, methods_group) %>% 
          summarise(mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE)), 
       aes(x = Method, y = mean, ymin = mean-sd, ymax = mean+sd, color = methods_group)) +
   geom_pointrange(show.legend = FALSE) +
   scale_color_viridis_d() +
   facet_wrap(N ~ ., ncol = 2, scales = "free_x") +
   coord_flip() +
   theme_bw() +
   labs(y = "Frobenius norm")



ggplot(combined_results_long %>% 
          dplyr::filter(Metric == "MAE_unique_values" & Method != "true_values") %>%  
          group_by(N, Method, Metric, methods_group) %>% 
          summarise(mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE)), 
       aes(x = Method, y = mean, ymin = mean-sd, ymax = mean+sd, color = methods_group)) +
   geom_pointrange(show.legend = FALSE) +
   scale_color_viridis_d() +
   facet_wrap(N ~ ., ncol = 2, scales = "free_x") +
   coord_flip() +
   theme_bw() +
   labs(y = "MAD unique values")



