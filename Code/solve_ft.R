if(!require(signal)){
  install.packages("signal")
  require(signal)
}

is_real_root <- function(all_roots){
  # User-defined rule for filtering real number.
  tor = 1e-7
  imag_num = Im(all_roots)
  index_logic = which(abs(imag_num) < tor)
  return(index_logic)
}

Mtx <- function(ts, x, e){
  num = size(ts)[1]
  if(num == 0){
    print('Test root is none.')
  }
  power_values = rep(0, e)
  for(i in 1:e)
    power_values[i] = 
  power_values = repmat(ts, 1, e)^(repmat(c(1:e), 1, num))
  Mts = M * t(power_values)
  results = Mts - repmat(x, 1, num)
  values = sum(results*results)
  return(values)
}

solve_ft <- function(M, x){
  e = size(M)[2]
  coefficients = array(0, dim = c(2*e, 1))
  
  # Calculate result_1(j, k) = M(:, j)'*M(:, k) for t^(j+k)
  result_1 = t(M) %*% M
  
  # Calculate result_2(k) = M(:, k)'*x for t^k
  result_2 = -2 * t(M) * x
  
  # Combine the final coefficients
  index_sum = repmat(c(1:e), e, 1) + repmat(matrix(c(1:e), nrow = length(c(1:e)), ncol = 1), 1, e)
  for(i in 2:(2*e)){
    coefficients[i] = coefficients[i] + sum(sum(result_1[index_sum == i]))
  }
  
  # Obtain the derivate coefficients and the corresponding roots
  derivate_coefficients_with_const = coefficients * matrix(c(1:(2*e)), nrow = length(c(1:(2*e))), ncol = 1)
  all_roots = roots(rev(derivate_coefficients_with_const))
  
  # Test obj. values of all real roots
  real_roots = all_roots[is_real_root(all_roots)]
  real_roots_objs = Mtx(real_roots, x, e)
  min_value = min(real_roots_objs)
  
  # There exists multiple minminizers.
  is_min = which(real_roots_objs - min_value == 0)
  minimizer_roots = real_roots[is_min]
  
  # Use the value-smallest minimizer as the calibration point.
  minimizer_root = min(minimizer_roots)
  if(nargout >= 2){
    minimum_value = min_value
  }
  
  result.final = data.frame(minimizer_root = minimizer_root, minimum_value = minimum_value)
  return(result.final)
}