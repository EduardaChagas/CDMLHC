if(!require(pracma)){
  install.packages("pracma")
  require(pracma)
}

gradient_regularizer <- function(tensor_M){
  grad = 2 * tensor_M
  return(grad)
}

obj_regularizer_R <- function(tensor_M){
  result = tensor_M * tensor_M
  value = sum(result)
  return(value)
}

batch_gradient_empirical_squared_dist <- function(squared_distance_values, labels){
  # Calculate the gradient of empirical loss w.r.t. squared_distance
  u = 2
  v = 8
  grad = labels * (2 * squared_distance_values[which(squared_distance_values > u)] * (squared_distance_values - u)) - 
        (1 - labels) * (2 * squared_distance_values[which(squared_distance_values < v)] * (v - squared_distance_values))
  return(grad)
}

obj_empirical_L <- function(tensor_M, X, X_hat, labels){
  d = size(tensor_M)[1]
  e = size(tensor_M)[2]
  m = size(tensor_M)[3]
  dd = size(X)[1]
  N = size(X)[2]
  empirical_loss_value = 0
  u = 2
  v = 8
  for(i in 1:N){
    empirical_loss_value = empirical_loss_value + labels[i] * max(0, squared_distance_value(tensor_M, X[,i], X_hat[,i]) - u)^2 +
                          (1- labels[i]) * max(0, v - squared_distance_value(tensor_M, X[, i], X_hat[, i]))^2
  }
  empirical_loss_value = empirical_loss_value / N
  return(empirical_loss_value)
}

cdml_training_sgd <- function(X, X_hat, labels, parameters){
  lambda = parameters$lambda
  c = parameters$c
  m = parameters$m
  batch_size = parameters$batch_size
  eta = parameters$eta
  tensor_M = parameters$tensor_M
  epoch = parameters$epoch
  
  d = size(X)[1]
  N = size(X)[2]
  obj_value = obj_empirical_L(tensor_M, X, X_hat, labels) + lambda * obj_regularizer_R(tensor_M)
  cat('Init: loss value = %f\n', obj_value)
  
  for(i in 1:epoch){
    indexs = randperm(N)
    batch_gap_check = 3
    batched = 0
    for(j in seq(from = 1, to = N, by = batch_size)){
      batch_indexs = indexs[c(j:min(j + batch_size - 1, N))]
      batch_num = length(batch_indexs)
      batch_grads = batch_gradient_squared_dist(tensor_M, X[, batch_indexs], X_hat[, batch_indexs])
      batch_grads_M = batch_grads[1]
      batch_squared_dists = batch_grads[2]
      batch_grads_squared_dist = batch_gradient_empirical_squared_dist(batch_squared_dists, labels[batch_indexs])
      temp = repmat(reshape(batch_grads_squared_dist, 1, 1, 1, batch_num), d, e, m, 1)
      grad_list = batch_grads_M * temp
      print(dim(grad_list))
      grad = sum(grad_list)
      tensor_M = tensor_M - eta*(grad + lambda*gradient_regularizer(tensor_M))
      batched  = batched + 1
      if(mod(batched, batch_gap_check)==0){
        obj_value = obj_empirical_L(tensor_M, X, X_hat, labels) + lambda * obj_regularizer_R(tensor_M)
        cat('epoch_', i, 'batch_', batched, ': loss value = ', obj_value, '\n')
      }
    }
  }
  return(tensor_M)
}