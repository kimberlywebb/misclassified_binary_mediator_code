sum_every_n <- function(x, n){
  vector_groups = split(x,
                        ceiling(seq_along(x) / n))
  sum_x = Reduce(`+`, vector_groups)

  return(sum_x)
}
