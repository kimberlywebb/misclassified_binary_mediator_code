sum_every_n1 <- function(x, n){
  vector_groups = split(x,
                        ceiling(seq_along(x) / n))
  sum_x = Reduce(`+`, vector_groups) + 1

  return(sum_x)
}
