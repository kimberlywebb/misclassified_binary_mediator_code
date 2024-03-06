# Simulation code:
## - Normal outcome (Y)
## - Medium misclassification rate (7% - 10%)
## - Predictive value weighting approach

################################################################################
# Load functions

#### UPDATE THIS FOR YOUR OWN COMPUTER! ####
function_directory <- "mediation_demo_for_github/functions/"
source_files <- list.files(function_directory, pattern = "*.R")
for (i in 1:length(source_files)) {
  source(paste0(function_directory, source_files[i]))
}

################################################################################
# Set save directory
save_directory <- "mediation_demo_for_github/results/"

################################################################################
# Load required packages
library(dplyr)
library(ggplot2)
library(stringr)
library(splines)

################################################################################
# Simulation settings

sample_size <- 10000
n <- sample_size

n_cat <- 2 # Number of categories in the binary mediator

# Data generation settings
x_mu <- 0
x_sigma <- 1
z_shape <- 1
z_scale <- 1
c_shape <- 1

# True parameter values (gamma terms set the misclassification rate)
true_beta <- matrix(c(1, -2, .5), ncol = 1)
true_gamma <- matrix(c(1.8, 1, -1.5, -1), nrow = 2, byrow = FALSE) 
true_theta <- matrix(c(1, 1.5, -2, -.2, .5), ncol = 1)

################################################################################
# Run simulations 

set.seed(123)

#### UPDATE THIS! ####
n_sim <- 10 # I chose a small number of simulations to make the run time short

# Initialize results vectors
my_results <- NULL
sim_indicator <- NULL

for(i in 1:n_sim){
  
  # Generate data
  data <- mediation_data_normalY_XM(sample_size,
                                    x_mu, x_sigma, z_shape, z_scale, c_shape,
                                    true_beta, true_gamma, true_theta)
  
  # Set starting values for the EM algorithm
  start_beta <- matrix(rep(1, 3), ncol = 1)
  start_gamma <- matrix(rep(1, 4), nrow = 2, ncol = 2)
  
  # Create matrix of true mediation model predictors
  mediation_model_predictors <- matrix(c(data[["x"]], data[["c"]]), ncol = 2,
                                       byrow = FALSE)
  
  # Run the COMBO EM algorithm for the true and observed mediation model
  mediation_results <- COMBO_EM_algorithm(data[["obs_mediator"]],
                                          mediation_model_predictors,
                                          data[["z"]],
                                          start_beta, start_gamma)
  
  # Save results
  predicted_beta <- matrix(mediation_results$Estimates[1:3], ncol = 1)
  predicted_gamma <- matrix(mediation_results$Estimates[4:7], 
                            ncol = 2, byrow = FALSE)
  
  # Create a matrix of observed mediator variables using dummy coding
  mstar_matrix <- matrix(c(ifelse(data[["obs_mediator"]] == 1, 1, 0),
                           ifelse(data[["obs_mediator"]] == 2, 1, 0)),
                         ncol = 2, byrow = FALSE)
  
  # Create matrix of predictors for the true mediator
  X_design <- matrix(c(rep(1, sample_size), data[["x"]], data[["c"]]),
                     ncol = 3, byrow = FALSE)
  
  # Generate probabilities for the true mediator value based on EM results
  pi_matrix <- pi_compute(predicted_beta, X_design, sample_size, n_cat)
  
  # Create matrix of predictors for the observed mediator
  Z_design <- matrix(c(rep(1, sample_size), data[["z"]]),
                     ncol = 2, byrow = FALSE)
  
  # Generate probabilities for observed mediator conditional on true mediator
  ## Based on EM results
  pistar_matrix <- pistar_compute(predicted_gamma, Z_design, sample_size, n_cat)
  
  # Estimate sensitivity and specificity
  sensitivity <- pistar_matrix[1:sample_size, 1]
  specificity <- pistar_matrix[(sample_size + 1):(2 * sample_size), 2]
  
  # Organize data for model predicting the observed mediator
  mstar_model_data <- data.frame(x = data[["x"]], c = data[["c"]], z = data[["z"]],
                                 y = data[["outcome"]],
                                 mstar = data[["obs_mediator"]])
  
  # Code observed mediator as 0/1 not 2/1
  mstar_model_data$mstar_01 <- ifelse(mstar_model_data$mstar == 1, 1, 0)
  
  # Fit spline model for observed mediator based on x, c, y, z
  mstar_model <- glm(mstar_01 ~ ns(y, 4) + ns(x, 4) + ns(c, 4) +
                       ns(x*c, 4) + ns(y*c, 4) + ns(z, 4),
                     data = mstar_model_data, family = "binomial")
  
  # Predict observed mediators
  predictions <- predict(mstar_model, type = "response")
  
  # Ensure no exact 0 or 1 values
  sensitivity[predictions >= sensitivity] <- predictions[predictions >= sensitivity] + 0.001
  specificity[predictions <= (1-specificity)] <- 1 - predictions[predictions <= (1 - specificity)] + 0.001
  
  # Compute NPV and PPV
  term1 <- (sensitivity - 1) * predictions * (1 / (sensitivity * (predictions - 1)))
  term2 <- (specificity - 1) * (predictions - 1) * (1 / (specificity * predictions))
  det <- 1/(term1*term2-1)
  ppv_calc <- det * (term2 - 1)
  npv_calc <- det * (term1 - 1)
  
  ppv <- unname(ppv_calc)
  npv <- unname(npv_calc)
  
  # Duplicate the dataset
  actual_dataset <- data.frame(x = data[["x"]], c = data[["c"]], z = data[["z"]],
                               y = data[["outcome"]], m = 0,
                               mstar_01 = mstar_model_data$mstar_01)
  
  duplicate_dataset <- data.frame(x = data[["x"]], c = data[["c"]], z = data[["z"]],
                                  y = data[["outcome"]], m = 1,
                                  mstar_01 = mstar_model_data$mstar_01)
  
  doubled_data <- rbind(actual_dataset, duplicate_dataset)
  
  # Apply NPV and PPV weights
  doubled_data$w <- 0
  doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 1] <- ppv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 1) - sample_size]
  doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 1] <- 1 - ppv[doubled_data$m == 0 & doubled_data$mstar_01 == 1]
  doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 0] <- 1 - npv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 0) - sample_size]
  doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 0] <- npv[doubled_data$m == 0 & doubled_data$mstar_01 == 0]
  
  # Fit weighted regression to estimate theta
  weighted_outcome_model <- lm(y ~ x + c + m + x*m, weights = w,
                               data = doubled_data)
  summary(weighted_outcome_model)
  
  # Save results and reoder theta
  results_i <- c(c(predicted_beta), c(predicted_gamma),
                 c(unname(coefficients(weighted_outcome_model))))
  
  my_results <- c(my_results, results_i)
  sim_indicator <- c(sim_indicator, rep(i, 12))
  
  print(paste0("Done with simulation ", i, "!"))
}

sim_results <- data.frame(estimate = my_results,
                          sim = rep(1:n_sim, each = 12),
                          parameter = rep(c("beta_0", "beta_x", "beta_c",
                                            "gamma_11", "gamma_21",
                                            "gamma_12", "gamma_22",
                                            "theta_0", "theta_x", "theta_c",
                                            "theta_m", "theta_xm"),
                                          n_sim),
                          true_value = rep(c(c(true_beta), c(true_gamma),
                                             c(true_theta)[c(1,2,4,3,5)]), n_sim))
# Summarize results
summarized_sim_results <- sim_results %>%
  mutate(bias = estimate - true_value) %>%
  group_by(parameter) %>%
  dplyr::summarise(mean_estimate = mean(estimate),
                   rmse_estimate = sqrt(((mean(estimate) - mean(true_value))^2) + (sd(estimate)^2)),
                   bias = mean(bias)) %>%
  ungroup()  %>%
  mutate(true_value = c(1, .5, -2, 1.8, -1.5, 1, -1, 1, -.2, -2, 1.5, .5))

summarized_sim_results

# Plot results
ggplot(data = sim_results) +
  geom_histogram(aes(x = estimate, fill = parameter),
                 color = "grey40", bins = 10) +
  geom_vline(aes(xintercept = true_value)) +
  theme_minimal() +
  facet_wrap(~parameter, scales = "free")

save(sim_results, summarized_sim_results,
     file = paste0(save_directory, "PVW_NormalY_MedMisclassification_results.RData"))
