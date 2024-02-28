# Function for response probabilities following a logistic model.
pi_compute <- function(true_beta, # Parameter values
                       X, # Design matrix
                       n, # Sample size
                       n_cat){
  
  exp_xb = exp(X %*% true_beta)
  pi_result = exp_xb[,1] / rep(sum_every_n1(exp_xb[,1], n), n_cat - 1)

  pi_matrix = matrix(c(pi_result, 1 - pi_result),
                     ncol = n_cat, nrow = n,
                     byrow = FALSE)
  return(pi_matrix)
}
