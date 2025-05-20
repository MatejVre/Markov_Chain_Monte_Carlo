library(ggplot2)
integrand <- function(x){
  (x ^ (-3/4)) * exp(-x)
} 

linspan = seq(0, 1, 0.01)

plot_df = data.frame(
  "X" = linspan,
  "Y" = integrand(linspan)
)

ggplot(data=plot_df, aes(x=X, y=Y)) +
  geom_line() +
  labs(
    x="x",
    y="y"
  ) +
  theme(
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave("integrand_plot.pdf", width = 6, height = 4)
#estimating the integral using uniform distribution
for(i in 1:10){
  n = 10^7
  prob = 1
  xs = runif(n, 0, 1)
  
  vals <- integrand(xs)
  mean <- 1/(n*prob) * sum(vals)
  std <- sd(vals)
  cat("Mean: ", mean, "+/-", std, "\n")
}
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

for (i in 1:10){
  xs_u <-  runif(n, 0, 1)
  xs_real <- inverse_cdf(xs_u)
  
  vals2 <- integrand(xs_real)/pdf(xs_real)
  
  eman <- mean(vals2)
  std <- sd(vals2)
  
  cat("Mean: ", mean, "+/-", std, "\n")
}

