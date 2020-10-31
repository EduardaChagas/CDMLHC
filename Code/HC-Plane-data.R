require(fftw)
source("Bandt-Pompe.R")

RNG.data <- function(){
  filePath <- read.table("/home/eduarda/Desktop/Repositories/ConfidenceRegions/Random_/Random_50k_D3-T1.dat", header=TRUE)
  HC <- matrix(data = filePath$x, nrow = 100, ncol = 3, byrow = TRUE, dimnames = NULL)
  HC <- data.frame(H = HC[,1], 
                   C = HC[,2], 
                   D = factor(rep(3, 100)))
  write.csv(HC, "White-noise.csv")
}

log.data <- function(){
  j = 1
  n = 50000
  n.series = 100
  x = seq(0, 2*pi, length.out = n*n.series)
  series.monotonic = log(x + 0.1)
  
  HC = data.frame(H = numeric(100), 
                  C = numeric(100), 
                  D = factor(rep(3, 100)))
  
  for(i in 1:n.series){
    time.series = series.monotonic[j:(j+n-1)]
    probs = bandt.pompe(time.series, 3, 1)
    HC$H[i] = shannon.entropy.normalized(probs)
    HC$C[i] = Ccomplexity(probs)
    j = j + n
  }
  write.csv(HC, "log-data.csv")
}

sin.data <- function(){
  j = 1
  n = 50000
  n.series = 100
  x = seq(0, 2*pi, length.out = n*n.series)
  series.periodic <- sin(2*x) * cos(2*x)
  
  HC = data.frame(H = numeric(100), 
                  C = numeric(100), 
                  D = factor(rep(3, 100)))
  
  for(i in 1:n.series){
    time.series = series.periodic[j:(j+n-1)]
    probs = bandt.pompe(time.series, 3, 1)
    HC$H[i] = shannon.entropy.normalized(probs)
    HC$C[i] = Ccomplexity(probs)
    j = j + n
  }
  write.csv(HC, "sin-data.csv")
}

series.generator.fk <- function(pp, y, n, k){
  Series <- vector(mode="numeric")
  filtro <- (1:n)^-(k/2)
  filtro <- filtro / sum(filtro)
  y1 <- y * filtro    
  x1 <- IFFT(y1, plan=pp)  
  Series <- c(Re(x1)) 
  Series
}

series.fk.data <- function(){
  n = 50000
  set.seed(seed = 1234567890, kind = "Mersenne-Twister")
  x = rnorm(n)
  x = x - mean(x)
  pp = planFFT(n)
  y = FFT(x, plan=pp)
  
  HC = data.frame(H = numeric(120), 
                  C = numeric(120),  
                  K = numeric(120), 
                  D = factor(rep(3, 120)))
  
  for(i in 1:120) {
    if(i <= 20){
      time.series = series.generator.fk(pp, y, n, .5)
      HC$K[i] = .5
    }
    else if(i <= 40){
      time.series = series.generator.fk(pp, y, n, 1)
      HC$K[i] = 1
    }
    else if(i <= 60){
      time.series = series.generator.fk(pp, y, n, 1.5)
      HC$K[i] = 1.5
    }
    else if(i <= 80){
      time.series = series.generator.fk(pp, y, n, 2)
      HC$K[i] = 2
    }
    else if(i <= 100){
      time.series = series.generator.fk(pp, y, n, 2.5)
      HC$K[i] = 2.5
    }
    else{
      time.series = series.generator.fk(pp, y, n, 3)
      HC$K[i] = 3
    }
    probs = bandt.pompe(time.series, 3, 1)
    HC$H[i] = shannon.entropy.normalized(probs)
    HC$C[i] = Ccomplexity(probs)
  }
  write.csv(HC, "series-fk-data.csv")
}
series.fk.data()
