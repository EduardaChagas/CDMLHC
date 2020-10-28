source('cdml_training_sgd.R')

cdml_hc <- function(){
  train.dataset = read.csv('../Data/train.csv', header=TRUE, sep=",")
  X = as.matrix(train.dataset[c('x_1', 'x_2', 'x_3', 'x_4')])
  X_hat = as.matrix(train.dataset[c('hat.x_1', 'hat.x_2', 'hat.x_3', 'hat.x_4')])
  labels = as.matrix(train.dataset['label'])
  parameters = data.frame( lambda = 0.001,
                           c = 2,
                           m = 2,
                           batch_size = 32,
                           eta = 0.01,
                           epoch = 100)
  tensor_M = array(1, dim = c(4, 2, 2))
  
  matrix.result = cdml_training_sgd(X, X_hat, labels, parameters, tensor_M)
  print(matrix.result)
}
