library(mcmcse)

# i will first mess around with markov chains to see if i understand it
unif_proposal_sampling <- function(x, delta) {
  runif(1, x - delta, x + delta)
}

unif_proposal_density <- function(x, x_condition, delta) {
  dunif(x, min=x_condition - delta, max=x_condition + delta)
}

my_dist <- function(x){return(dnorm(x, 0, 1))}

metropolis_hastings <- function(x, p, proposal_sampler, proposal_density, lower_bound= -Inf, upper_bound = Inf, ...){
  x_new <- proposal_sampler(x, ...)
  
  if (any(x_new <= lower_bound) || any(x_new >= upper_bound)) {
    return(x)  # reject proposal immediately
  }
  #correction step
  #alpha <- min((p(x_new) / p(x)), 1) #this is Metropolis only because the distribution is the same
  #alpha <- min((p(x_new) * proposal_density(x_new, x, ...) / p(x) * proposal_density(x, x_new, ...)), 1) #Metropolis-Hastings update step
  #it all makes sense now, Metropolis-Hastings with a symmetric proposal is just a metropolis algorithm
  
  log_alpha <- log(p(x_new)) + log(proposal_density(x, x_new, ...)) - log(p(x)) - log(proposal_density(x_new, x, ...))
  
  if (!is.finite(log_alpha)){
    alpha <- 0
  }else{
    alpha <- min(exp(log_alpha), 1)
  }
  
  u <- runif(1)
  
  if (u > alpha){x_new <- x}
 # if (x_new > upper_bound | x_new < lower_bound) { x_new <- x } #left the domain
  
  return(x_new)
}

#m <- 20000
#d <- 0.78
#samp <- c(1)

#for (i in 1:m){
#  x_new <- metropolis_hastings(samp[length(samp)], my_dist, unif_proposal_sampling, unif_proposal_density, delta=d)
  
#  samp <- c(samp, x_new)
#}

#hist(samp)
#mean(samp)
#sd(samp) / sqrt(m)
#acf(samp)
#ess(samp)
#plot(samp, type="l")

seen_data <- c(0.3, 0.5, 0.75, 0.4)

weibull_likelihood <- function(data, alpha, eta){
  return(alpha*eta*data^(alpha - 1)*exp(-(data^(alpha))*eta))
}

prior_density <- function(alpha, eta){
  return(exp(-alpha-2*eta)*eta)
}

posterior <- function(params, data=seen_data){
  alpha <- params[1]
  eta <- params[2]
  likelihood <- prod(weibull_likelihood(data, alpha, eta))
  prior <- prior_density(alpha, eta)
  return(likelihood*prior)
}

library(mvtnorm)
mnorm_proposal_sampling <- function(x, sigma){
  rmvnorm(1, x, sigma)
}

mnorm_proposal_density <- function(x, x_condition, sigma){
  dmvnorm(x, x_condition, sigma)
}

#M-H using multivariate normal distribution------------------------------------------------------------------
set.seed(42)
sigma_matrix <- matrix(c(1, 0, 0, 1), ncol = 2)
chains = list()
for (j in 1:5){
  m <- 1000
  samples <- data.frame(X = 0.5, Y = 0.5) # starting value
  for (i in 1:m) {
    x <- samples[nrow(samples),]
    c(x$X, x$Y)
    x_new <- metropolis_hastings(c(x$X, x$Y), posterior, mnorm_proposal_sampling, mnorm_proposal_density, lower_bound = 0, upper_bound = Inf, sigma = sigma_matrix)
    samples <- rbind(samples, data.frame(X = x_new[1], Y = x_new[2]))
  }
  chains[[j]] <- samples
}

#It would appear that setting this is extremely difficult, i have been searching the
#sigma matrix parameter space for the better part of an hour. *angry_hello_kitty_emoji*
plot_autocorrelation <- function(chains, filename){
  
  #windows(width = 40, height = 20)
  grDevices::pdf(filename, width = 16, height = 10)
  par(mfrow=c(5,4))
  
  for (i in 1:5){
    samples <- chains[[i]]
    alpha_samples <- samples$X
    eta_samples <- samples$Y
    
    acf(alpha_samples, main = paste("Alpha Chain", i, "- Autocorrelation"))
    plot(alpha_samples, type="l", main = paste("Alpha Chain", i, "- Traceplot"))
    acf(eta_samples, main = paste("Eta Chain", i, "- Autocorrelation"))
    plot(eta_samples, type="l", main = paste("Eta Chain", i, "- Traceplot"))
  }
  
  dev.off()
}

plot_autocorrelation2 <- function(chains, filename){
  
  #windows(width = 40, height = 20)
  grDevices::pdf(filename, width = 16, height = 16)
  par(mfrow=c(5,2), cex.main = 2.0, cex.lab = 1.7, cex.axis = 1.7)
  
  for (i in 1:5){
    samples <- chains[[i]]
    alpha_samples <- samples$X
    eta_samples <- samples$Y
    
    acf(alpha_samples, main = paste("Alpha Chain", i, "- Autocorrelation"))
    acf(eta_samples, main = paste("Eta Chain", i, "- Autocorrelation"))
  }
  
  dev.off()
}

plot_trace <- function(chains, filename){
  
  #windows(width = 40, height = 20)
  grDevices::pdf(filename, width = 16, height = 16)
  par(mfrow=c(5,2), cex.main = 2.0, cex.lab = 1.7, cex.axis = 1.7)
  
  for (i in 1:5){
    samples <- chains[[i]]
    alpha_samples <- samples$X
    eta_samples <- samples$Y
    
    plot(alpha_samples, type="l", main = paste("Alpha Chain", i, "- Traceplot"))
    plot(eta_samples, type="l", main = paste("Eta Chain", i, "- Traceplot"))
  }
  
  dev.off()
}

ess_for_all <- function(chains){
  ess_alphas <- c()
  ess_etas <- c()
  for (ch in chains){
    alphas <- ch$X
    etas <- ch$Y
    ea <- ess(alphas)
    ee <- ess(etas)
    ess_alphas <- c(ess_alphas, ea)
    ess_etas <- c(ess_etas, ee)
    cat("Alpha ess: ", ea, "|", "Eta ess: ", ee, "\n")
  }
  cat("Mean alpha ess: ", mean(ess_alphas), "|", "Mean eta ess: ", mean(ess_etas), "\n")
}

acceptance_rate_for_all <- function(chains){
  for (ch in chains){
    alphas <- ch$X
    etas <- ch$Y
    cat("Alpha rejection rate: ", sum(duplicated(alphas))/m, "|", "Eta rejection rate: ", sum(duplicated(etas)) / m, "\n")
  }
}

chain_mean_and_variance <- function(chains){

  for (ch in chains){
    alphas <- ch$X
    etas <- ch$Y
    cat("Mean Alpha: ", mean(alphas), "+/-", mcse(alphas)$se^2, "|", "Mean Eta: ", mean(etas), "+/-", mcse(etas)$se^2, "\n")
  }
}

joined_chain_mean_and_variance <-  function(chains){
  alphas <- c()
  etas <- c()
  for (ch in chains){
    alphas <- c(alphas, ch$X)
    etas <- c(etas, ch$Y)
  }
  cat("Mean Alpha: ", mean(alphas), "+/-", var(alphas), "|", "Mean Eta: ", mean(etas), "+/-", var(etas), "\n")
}

chain_cdf <- function(x_alpha=1.3, x_eta=1.3, chains=chains){
  for (ch in chains){
    alphas <- ch$X
    etas <- ch$Y
    
    p_alpha <- length(alphas[alphas >= x_alpha]) / length(alphas)
    p_eta <- length(etas[etas >= x_eta]) / length(etas)
    
    cat("Probability: ", p_alpha * p_eta, "\n")
  }
}

joined_chain_cdf <- function(x_alpha=1.3, x_eta=1.3, chains=chains){
  alphas <- c()
  etas <- c()
  for (ch in chains){
    alphas <- c(alphas, ch$X)
    etas <- c(etas, ch$Y)
    
  }
  p_alpha <- length(alphas[alphas >= x_alpha]) / length(alphas)
  p_eta <- length(etas[etas >= x_eta]) / length(etas)
  
  cat("Probability: ", p_alpha * p_eta, "\n")
}
#Diagnostics and results------------------------------------------------------------------
plot_autocorrelation2(chains, "multi_autocorrelation.pdf")
plot_trace(chains, "multi_trace.pdf")
ess_for_all(chains)
acceptance_rate_for_all(chains)
chain_mean_and_variance(chains)
joined_chain_mean_and_variance(chains)
chain_cdf(chains=chains)
joined_chain_cdf(chains=chains)


alpha_seq <- seq(0.1, 5, length.out = 100)
eta_seq <- seq(0.1, 5, length.out = 100)

# Compute posterior values over the grid
posterior_grid <- outer(alpha_seq, eta_seq, Vectorize(function(a, e) posterior(c(a, e))))

par(mfrow = c(1, 1))
contour(alpha_seq, eta_seq, posterior_grid,
        xlab = "Alpha", ylab = "Eta", main = "Posterior Contour with Chain")
samples <- chains[[1]]

lines(samples$X, samples$Y, col = "red", lwd = 2)      # chain path
points(samples$X, samples$Y, col = "blue", pch = 20)


#M-H using the other distribution########################

exp_proposal_sampling <- function(params){
  alpha <- params[1]
  eta <- params[2]
  
  alpha_prime <- rexp(1, 1/alpha)
  eta_prime <- rexp(1, 1/eta)
  
  return(c(alpha_prime, eta_prime))
}

exp_proposal_density <- function(params, params_condition){
  alpha <- params[1]
  eta <- params[2]
  alpha_condition <- params_condition[1]
  eta_condition <- params_condition[2]

  return(dexp(alpha, 1/alpha_condition) * dexp(eta, 1/eta_condition))  
}

set.seed(42)
sigma_matrix <- matrix(c(0.2, 0.11, 0.11, 0.18), ncol = 2)
chains_other = list()
for (j in 1:5){
  m <- 1000
  samples <- data.frame(X = 0.5, Y = 0.5) # starting value
  for (i in 1:m) {
    x <- samples[nrow(samples),]
    c(x$X, x$Y)
    x_new <- metropolis_hastings(c(x$X, x$Y), posterior, exp_proposal_sampling, exp_proposal_density, lower_bound = 0, upper_bound = Inf)
    samples <- rbind(samples, data.frame(X = x_new[1], Y = x_new[2]))
  }
  chains_other[[j]] <- samples
}

plot_autocorrelation(chains_other)
ess_for_all(chains_other)
acceptance_rate_for_all(chains_other)
chain_mean_and_variance(chains_other)
chain_cdf(chains=chains_other)
