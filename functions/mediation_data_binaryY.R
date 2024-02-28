mediation_data_binaryY <- function(sample_size,
                                   x_mu, # Mean of X ~ Normal
                                   x_sigma, # SD of X ~ Normal
                                   z_shape, # Shape for Z ~ Gamma
                                   z_scale, # Scale for Z ~ Gamma
                                   c_shape, # Shape for C ~ Gamma
                                   true_beta, true_gamma, true_theta){
  
  n_cat <- 2 # Number of categories in mediator
  
  # Generate X
  x <- rnorm(sample_size, x_mu, x_sigma)
  x_matrix <- matrix(c(rep(1, sample_size),
                       x),
                     nrow = sample_size, byrow = FALSE)

  # Generate Z
  z <- rgamma(sample_size, z_shape)
  z_matrix <- matrix(c(rep(1, sample_size),
                       z),
                     nrow = sample_size, byrow = FALSE)

  # Generate C
  c <- rgamma(sample_size, c_shape)
  c_matrix <- matrix(c(rep(1, sample_size),
                       c),
                     nrow = sample_size, byrow = FALSE)
  
  # Create matrix of predictors for the true mediator
  predictor_matrix <- cbind(x_matrix, c_matrix[,2])
  
  # Generate probabilities for the true mediator value
  pi_matrix <- pi_compute(true_beta, predictor_matrix, sample_size, n_cat)

  # Generate true mediator variable based on probabilities in pi_matrix
  true_M <- rep(NA, sample_size)
  for(i in 1:sample_size){
    true_M[i] = which(rmultinom(1, 1, pi_matrix[i,]) == 1)
  }
  
  # Generate probabilities for observed mediator conditional on true mediator
  pistar_matrix <- pistar_compute(true_gamma, z_matrix, sample_size, n_cat)
  
  # Generate observed mediator variable based on conditional probabilities in pistar_matrix
  obs_Mstar <- rep(NA, sample_size)
  for(i in 1:sample_size){
    true_j = true_M[i]
    obs_Mstar[i] = which(rmultinom(1, 1,
                                   pistar_matrix[c(i,sample_size + i),
                                                 true_j]) == 1)
  }
  
  # Create a matrix of observed mediator variables using dummy coding
  obs_Mstar_reps <- matrix(rep(obs_Mstar, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix <- matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                            byrow = FALSE)
  obs_Mstar_matrix <- 1 * (obs_Mstar_reps == category_matrix)

  # Generate probabilities for the outcome
  pitilde_matrix <- prob_yi(true_theta, true_M, x, c, sample_size, n_cat)

  # Generate binary outcome based on probabilities in pitilde_matrix
  outcome <- rep(NA, sample_size)
  for(i in 1:sample_size){
    outcome[i] = which(rmultinom(1, 1, pitilde_matrix[i,]) == 1)
  }
  
  # Convert outcome to 0/1 instead of 2/1 variable coding
  outcome_01 = ifelse(outcome == 1, 1, 0)

  # Organize data for output
  data_output <- list(obs_mediator = obs_Mstar,
                      outcome = outcome_01,
                      true_mediator = true_M,
                      x = x,
                      z = z,
                      c = c,
                      x_design_matrix = x_matrix,
                      z_design_matrix = z_matrix,
                      c_design_matrix = c_matrix)

  return(data_output)

}

