source('arc_length.R')
source('solve_ft.R')

squared_distance_value <- function(tensor_M, x, x_hat){
  d = size(tensor_M)[1]
  e = size(tensor_M)[2]
  m = size(tensor_M)[3]
  value = 0
  for(i in 1:m){
    value = value + 
            arc_length(tensor_M[,,i], 0, 1)^2 *
            arc_length(tensor_M[,,i], solve_ft(tensor_M[,,i], x), solve_ft(tensor_M[,,i], x_hat))^2
            
  }
  return(value)
}