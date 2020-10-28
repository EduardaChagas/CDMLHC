if(!require(pracma)){
  install.packages("pracma")
  require(pracma)
}

if(!require(matconv)){
  install.packages("matconv")
  require(matconv)
}

arc_length <- function(M, A, B){
  # The number of intervals
  L = 1000
  
  # Calculate arc length by interval approximation
  d = dim(M)[1]
  e = dim(M)[2]
  a = min(A, B)
  b = max(A, B)
  delta_t = (b - a)/L
  inter_values = a + c(0:(L-1)) * delta_t
  
  # Pre-calcualte power values at intervals
  power_inter_values = repmat(inter_values, e, 1)^(repmat(matrix(c(0:(e-1)), nrow = length(c(0:(e-1))), ncol = 1), 1, length(inter_values)))
  
  # Convert poly. computations to linear projections
  derivate_coefficients_with_const = M * repmat(c(1:e), d, 1)
  result = derivate_coefficients_with_const %*% power_inter_values
  value = sum(sqrt(colSums(result * result))) * delta_t
  return(value)
}
