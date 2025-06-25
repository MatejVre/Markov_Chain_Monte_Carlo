library(mc2d)

a = 0
b = 10
c = 100

Ex = (a + 4*b + c)/6

pert <- function(x, a=0, b=10, c=100){
  alpha = 1 + 4*((b-a)/(c-a))
  bet = 1 + 4*((c-b)/(c-a))
  
  numerator = ((x - a) ^ (alpha-1)) * ((c-x)^(bet-1))
  denominator = beta(alpha, bet)*((c-a)^(alpha+bet-1))
  
  return(numerator/denominator)
}

trapezoidal_rule_for_expected_value <- function(lower_bound, upper_bound, fun, N){
  
  cum_sum <- 0
  h = (upper_bound - lower_bound) / N
  
  for (x in seq(lower_bound, upper_bound, h)){
    fx <- x * fun(x)
    if (x == lower_bound | x == upper_bound){
      cum_sum <- cum_sum + fx
    }
    else{
      cum_sum <- cum_sum + 2*fx
    }
  }
  return((h / 2)*cum_sum)
}

#It takes 271 evaluations to get the value approximated to 4 decimal places
res <- trapezoidal_rule_for_expected_value(a, c, pert, 271)
cat(sprintf("approx = %.8f", res))
print(abs(res-Ex)<0.00005)


#trying the CLT
n = 1e8
xs = runif(n, a, c)
prob = 1/c-a

f_samples <- xs * pert(xs) * (c - a)

ex_est <- mean(f_samples)

#ex2_est <- mean(f_samples * (c - a))

sample_sd <- sd(f_samples)

required_n <- ceiling((1.96 * sample_sd / 0.005)^2)

#We require 55939373 samples to estimate E[X] to 2 decimal places 95% of the time.

#Doing the numerical examples but might switch to Python if it'll be too slow.

for(i in 1:200){
  
}
