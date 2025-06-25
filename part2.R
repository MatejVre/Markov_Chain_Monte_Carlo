library(ggplot2)
integrand <- function(x){
  (x ^ (-3/4)) * exp(-x)
} 

linspan = seq(0, 1, 0.01)

plot_df = data.frame(
  "X" = linspan,
  "Y" = integrand(linspan),
  "Y2" = pdf(linspan)
)

ggplot(data = plot_df) +
  geom_line(aes(x = X, y = Y, color = "Integrand")) +
  geom_line(aes(x = X, y = Y2, color = "Proposal PDF")) +
  labs(
    x = "x",
    y = "y",
    color = "Legend"
  ) +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

ggsave("integrand_plot.pdf", width = 6, height = 4)
#estimating the integral using uniform distribution
uniform_means <- c()
for(i in 1:10){
  n = 10^7
  prob = 1
  xs = runif(n, 0, 1)
  
  vals <- integrand(xs)
  mean <- 1/(n*prob) * sum(vals)
  uniform_means <- c(uniform_means, mean)
}
cat("Mean: ", mean(uniform_means), "+/-", sd(uniform_means), "\n")
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

better_means <- c()
for (i in 1:10){
  xs_u <-  runif(n, 0, 1)
  xs_real <- inverse_cdf(xs_u)
  
  vals2 <- integrand(xs_real)/pdf(xs_real)
  
  better_means <- c(better_means, mean(vals2))
  
  
}
cat("Mean: ", mean(better_means), "+/-", sd(better_means), "\n")
