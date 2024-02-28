COMBO_EM_function <- function(param_current,
                              obs_Y_matrix, X, Z,
                              sample_size, n_cat){

  beta_current = matrix(param_current[1:ncol(X)], ncol = 1)
  gamma_current = matrix(c(param_current)[(ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))],
                         ncol = n_cat, byrow = FALSE)

  probabilities = pi_compute(beta_current, X, sample_size, n_cat)
  conditional_probabilities = pistar_compute(gamma_current, Z, sample_size, n_cat)

  weights = COMBO_weight(ystar_matrix = obs_Y_matrix,
                         pistar_matrix = conditional_probabilities,
                         pi_matrix = probabilities,
                         sample_size = sample_size, n_cat = n_cat)

  Ystar01 = obs_Y_matrix[,1]
  fit.gamma1 <- suppressWarnings( stats::glm(Ystar01 ~ . + 0, as.data.frame(Z),
                           weights = weights[,1],
                           family = "binomial"(link = "logit")) )
  gamma1_new <- unname(coefficients(fit.gamma1))

  fit.gamma2 <- suppressWarnings( stats::glm(Ystar01 ~ . + 0, as.data.frame(Z),
                           weights = weights[,2],
                           family = "binomial"(link = "logit")) )
  gamma2_new <- unname(coefficients(fit.gamma2))

  fit.beta <- suppressWarnings( stats::glm(weights[,1] ~ . + 0, as.data.frame(X),
                         family = stats::binomial()) )
  beta_new <- unname(coefficients(fit.beta))

  param_new = c(beta_new, gamma1_new, gamma2_new)
  param_current = param_new
  return(param_new)

}
