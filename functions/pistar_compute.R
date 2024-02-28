# Function for response probabilities following a logistic model with a
## conditional outcome
pistar_compute <- function(gamma, # Parameter values
                           Z, # Design matrix
                           n, # Sample size
                           n_cat){

  exp_zg = exp(Z %*% gamma)
  pi_denominator = apply(exp_zg, FUN = sum_every_n1, n, MARGIN = 2)
  pi_result = exp_zg / rbind(pi_denominator)

  pistar_matrix = rbind(pi_result,
                        1 - apply(pi_result,
                                  FUN = sum_every_n, n = n,
                                  MARGIN = 2))

  return(pistar_matrix)
}
