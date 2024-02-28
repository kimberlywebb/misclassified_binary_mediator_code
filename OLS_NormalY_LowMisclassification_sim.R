# Simulation code:
## - Normal outcome (Y)
## - Low misclassification rate (2% - 5%)
## - Ordinary least squares correction approach

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
true_gamma <- matrix(c(3, 2, -2, -2.5), nrow = 2, byrow = FALSE)
true_theta <- matrix(c(1, 1.5, -2, -.2), ncol = 1)
################################################################################
# Run simulations 

set.seed(123)

#### UPDATE THIS! ####
n_sim <- 10 # I chose a small number of simulations to make the run time short

# Initialize results vectors
my_results <- NULL
sim_indicator <- NULL

for(i in 1:n_sim){
  
  data <- mediation_data_normalY(sample_size,
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
  
  # Compute the observed mediator prevalence
  prevalence <- length(which(data[["obs_mediator"]] == 1)) / sample_size
  
  # Compute average misclassification rates
  pistar12 <- pistar_matrix[1:sample_size, 2]
  pistar21 <- pistar_matrix[(sample_size + 1):(2 * sample_size), 1]
  
  # Compute correction parameters from Nguimkeu, Rosenman, and Tennekoon (2021)
  theta_Nguimkeu <- (pistar12 + pistar21) / (1 - pistar12 - pistar21)
  squiggle_Nguimkeu <- 1 - (((prevalence - pistar12)*(1 - pistar21 - prevalence)) / 
                              ((1 - pistar12 - pistar21)*(1 - prevalence)*prevalence))
  
  # Compute covariances for the correction
  m_matrix <- matrix(ifelse(data[["obs_mediator"]] == 1, 1, 0), ncol = 1)
  sd_dd <- cov(m_matrix)
  
  predictor_matrix <- matrix(c(data[["x"]], data[["c"]]), ncol = 2, byrow = FALSE)
  sd_xx <- cov(predictor_matrix)
  
  sd_xd <- cov(predictor_matrix, m_matrix)
  sd_dx <- cov(m_matrix, predictor_matrix)
  
  y_matrix <- matrix(data[["outcome"]], ncol = 1)
  sd_yd <- cov(y_matrix, m_matrix)
  sd_yx <- cov(y_matrix, predictor_matrix)
  
  block1_dd <- (1 - median(squiggle_Nguimkeu)) * sd_dd[1,1]
  block1_xd <- (1 + median(theta_Nguimkeu)) * sd_xd
  block_1_matrix <- matrix(c(block1_dd, block1_xd[,1],
                             sd_dx[,1], sd_xx[,1],
                             sd_dx[,2], sd_xx[,2]), byrow = FALSE,
                           nrow = 3)
  
  block_2_matrix <- matrix(c(sd_yd, sd_yx[1,]), ncol = 1)
  
  # Solve for the corrected parameters
  solve_param <- solve(block_1_matrix) %*% block_2_matrix
  solve_param
  
  # Compute the intercept
  intercept <- mean(data[["outcome"]]) -
    ((solve_param[1,1]) * (mean(m_matrix[,1]) - median(pistar12)) / (1 - median(pistar12) - median(pistar21))) -
    t(colMeans(predictor_matrix) %*% solve_param[-1,])
  
  # Intercept not estimated in "solve_param".
  results_i <- c(c(predicted_beta), c(predicted_gamma),
                 intercept, solve_param[c(2, 1, 3), 1])
  
  my_results <- c(my_results, results_i)
  sim_indicator <- c(sim_indicator, rep(i, 11))
  
  print(paste("Done with simulation", i))
}

sim_results <- data.frame(estimate = na.omit(my_results),
                          sim = rep(1:n_sim, each = 11),
                          parameter = rep(c("beta_0", "beta_x", "beta_c",
                                            "gamma_11", "gamma_21",
                                            "gamma_12", "gamma_22",
                                            "theta_0", "theta_x", "theta_m",
                                            "theta_c"),
                                          n_sim),
                          true_value = rep(c(c(true_beta), c(true_gamma),
                                             c(true_theta)), n_sim))

# Summarize results
summarized_sim_results <- sim_results %>%
  mutate(bias = estimate - true_value) %>%
  group_by(parameter) %>%
  dplyr::summarise(mean_estimate = mean(estimate),
                   rmse_estimate = sqrt(((mean(estimate) - mean(true_value))^2) + (sd(estimate)^2)),
                   bias = mean(bias)) %>%
  ungroup()  %>%
  mutate(true_value = c(1, .5, -2, 3, -2, 2, -2.5, 1, -.2, -2, 1.5))

summarized_sim_results

# Plot results
ggplot(data = sim_results) +
  geom_histogram(aes(x = estimate, fill = parameter),
                 color = "grey40", bins = 10) +
  geom_vline(aes(xintercept = true_value)) +
  theme_minimal() +
  facet_wrap(~parameter, scales = "free")

save(sim_results, summarized_sim_results,
     file = paste0(save_directory,
                   "OLS_NormalY_LowMisclassification_results.RData"))
