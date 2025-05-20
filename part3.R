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
  
  #correction step
  alpha <- min((p(x_new) / p(x)), 1) #this is metropolis only because the distribution is the same
  alpha <- min((p(x_new) * proposal_density(x_new, x, ...) / p(x) * proposal_density(x, x_new, ...)), 1) #metropolis hastings update step
  
  u <- runif(1)
  
  if (u > alpha){x_new <- x}
  
  if (x_new > upper_bound | x_new < lower_bound) { x_new <- x } #left the domain
  
  return(x_new)
}

m <- 10000
samp <- c(1)

for (i in 1:m){
  x_new <- metropolis_hastings(samp[length(samp)], my_dist, unif_proposal_sampling, unif_proposal_density, delta=0.5)
  
  samp <- c(samp, x_new)
}

hist(samp)
