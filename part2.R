integrand <- function(x){
  (x ^ (-3/4)) * exp(-x)
} 

linspan = seq(0, 1, 0.01)

plot(linspan, integrand(linspan))

#estimating the integral using uniform distribution

n = 10^7
prob = 1
xs = runif(n, 0, 1)

1/(n*prob) * sum(integrand(xs))

#estimating the integral using inversion sampling
pdf <- function(x){
  return(1/4*x^(-3/4))
}

cdf <- function(x){
  return(x^(1/4))
}

inverse_cdf <- function(x){
  return(x^4)
}

xs_u <-  runif(n, 0, 1)
xs_real <- inverse_cdf(xs_u)

1/n * sum(integrand(xs_real)/pdf(xs_real))


