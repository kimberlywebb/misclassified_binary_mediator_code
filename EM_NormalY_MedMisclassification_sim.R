# Simulation code:
## - Normal outcome (Y)
## - Medium misclassification rate (7% - 10%)
## - EM algorithm

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
  start_theta <- matrix(rep(1, 5), ncol = 1)
  start_sigma <- 1
  
  # Run the EM algorithm
  mediation_results <- mediation_EM_normalY_XM(Mstar = data[["obs_mediator"]],
                                               outcome = data[["outcome"]],
                                               x_matrix = data[["x"]],
                                               z_matrix = data[["z"]],
                                               c_matrix = data[["c"]],
                                               beta_start = start_beta,
                                               gamma_start = start_gamma,
                                               theta_start = start_theta,
                                               sigma_start = start_sigma,
                                               tolerance = 1e-7,
                                               max_em_iterations = 1500,
                                               em_method = "squarem")
  
  # Update results vector
  results_i <- mediation_results$Estimates
  my_results <- c(my_results, results_i)
  
  # Update simulation indicator vector
  sim_indicator <- c(sim_indicator, rep(i, 13))
  
  # Print simulation number
  print(paste0("Done with simulation ", i, "!"))
}

# Put results into data frame
sim_results <- data.frame(estimate = my_results,
                          sim = sim_indicator,
                          parameter = rep(c("beta_0", "beta_c", "beta_x",
                                            "gamma_11", "gamma_21",
                                            "gamma_12", "gamma_22",
                                            "theta_0", "theta_x", "theta_m",
                                            "theta_c", "theta_xm",
                                            "sigma"),
                                          n_sim),
                          true_value = rep(c(c(true_beta), c(true_gamma),
                                             c(true_theta), 1),
                                           n_sim))

# Summarise results
summarized_sim_results <- sim_results %>%
  mutate(bias = estimate - true_value) %>%
  group_by(parameter) %>%
  dplyr::summarise(mean_estimate = mean(estimate),
                   rmse_estimate = sqrt(((mean(estimate) - mean(true_value))^2)
                                        + (sd(estimate)^2)),
                   bias = mean(bias)) %>%
  ungroup()  %>%
  mutate(true_value = c(1, -2, .5, 1.8, -1.5, 1, -1, 1, 1, -.2, -2, 1.5, .5))

summarized_sim_results

# Plot results
ggplot(data = sim_results) +
  geom_histogram(aes(x = estimate, fill = parameter),
                 color = "grey40", bins = 10) +
  geom_vline(aes(xintercept = true_value)) +
  theme_minimal() +
  facet_wrap(~parameter, scales = "free")

# Save results
save(sim_results, summarized_sim_results,
     file = paste0(save_directory, "EM_NormalY_MedMisclassification_results.RData"))
