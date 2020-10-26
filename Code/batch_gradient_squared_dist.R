batch_gradient_squared_dist <- function(tensor_M, X, X_hat){
  # Caculate the batch gradients and squared distances.
  d = size(tensor_M)[1]
  e = size(tensor_M)[2]
  m = size(tensor_M)[3]
  
  dd = size(X)[1]
  batch_num = size(X)[2]
  grads = array(0, dim = c(d, e, m, batch_num))
  squared_dists = array(0, dim = c(batch_num, 1))
  
  foreach(i = 1:batch_num) %do%{
    for(j in 1:m){
      Mj = tensor_M[, , j]
      x = X[, i]
      x_hat = X_hat[, i]
      cali_x = solve_ft(Mj, x)
      cali_x_hat = solve_ft(Mj, x_hat)
      length_unit = arc_length(Mj, 0, 1)
      length_cali = arc_length(Mj, cali_x, cali_x_hat)
      grads[, , j, i] = 2 * length_unit * gradient_arc_length(Mj, 0, 1) * 
                        (length_cali^2) + (length_unit^2) * gradient_cali_squared_arc_length(Mj, x, x_hat)
      squared_dists[i] = squared_dists[i] + (length_unit^2) * (length_cali^2)
    }
  }
  result = data.frame(grads = grads, squared_dists = squared_dists)
  return(result)
}

gradient_arc_length <-function(Mi, A, B){
  d = size(Mi)[1]
  e = size(Mi)[2]
  L = 1000
  a = min(A, B)
  b = max(A, B)
  delta_t = (b - a)/L
  inter_values = a + c(0:(L-1))*delta_t
  
  # Pre-calcualte power values at intervals
  power_inter_values = repmat(inter_values, c, 1)^(repmat(t(c(0:(e-1))), 1, length(inter_values)))
  
  # Derivate
  derivate_coefficients_with_const = Mi * repmat(c(1:e), d, 1)
  result = derivate_coefficients_with_const * power_inter_values[c(1:e), ]
  result_1 = 1 / sqrt(sum(result * result))
  
  # Sum operation
  grad = delta_t * (result * diag(result_1) * t(power_inter_values))
  return(grad)
}

gradient_cali_squared_arc_length <- function(Mi, x, x_hat){
  d = size(Mi)[1]
  e = size(Mi)[2]
  mu = 1e-2
  delta_M = randn(d, e)
  mu_delta_M = mu * delta_M
  norm_delta = norm(delta_M)
  Mi_move = Mi + mu_delta_M
  
  cali_x = solve_ft(Mi, x)
  cali_x_hat = solve_ft(Mi, x_hat)
  cali_x_move = solve_ft(Mi_move, x)
  cali_x_hat_move = solve_ft(Mi_move, x_hat)
  grad = (delta_M/(mu*(norm_delta^2))) * 
         (arc_length(Mi_move, cali_x_move, cali_x_hat_move)^2 - arc_length(Mi, cali_x, cali_x_hat)^2)
  return(grad)
}