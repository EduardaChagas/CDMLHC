
if(!require(matconv)){
  install.packages("matconv")
  require(matconv)
}

matIn <- c("function value = arc_length(M, A, B)", 
           " L = 1000;",
           " [d, c] = size(M);",
           " a = min(A, B);",
           " b = max(A, B);",
           " delta_t = (b - a)/L;",
           " inter_values = a + (0:L-1)*delta_t;",
           " power_inter_values = repmat(inter_values, c, 1) ...
                  .^ repmat((0:c-1)', 1, length(inter_values));",
           " derivate_coefficients_with_const = M .* repmat(1:c, d, 1);",
           " result = derivate_coefficients_with_const * power_inter_values;",
           " value = sum(sqrt(sum(result .* result, 1))) * delta_t;",
           "end"
           )

mat2r(matIn, verbose = 0)$rCode
