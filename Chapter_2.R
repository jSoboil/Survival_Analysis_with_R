# ====================================================================================
# Basic Principles of Survival Analysis -----------------------------------
# ====================================================================================
library(survival)

# In R the functions “dweibull” and “pweibull” compute the p.d.f. and c.d.f., respectively, 
# of the Weibull distribution. These functions use the arguments “shape” and “scale” to 
# represent the parameters alpha and 1/lambda, respectively. To obtain the survival function, 
# we can specify “lower.tail = F” as an option in the “pweibull” function. For example, we 
# can plot the Weibull survival function with alpha = 1 and lambda = .03 by first defining a
# function “weibSurv” with these parameters and then using the “curve” function to plot the 
# curve as follows
weibSurv <- function(t, shape, scale) {
 pweibull(t, shape = shape, scale = scale, lower.tail = F)
}

curve(weibSurv(x, shape = 1, scale = 1 / .03), from = 0, to = 80,
      ylim = c(0, .08), ylab = "Survival probability", xlab = "Time", col = "darkblue")

# To plot the hazard function with this shape and scale, we can use the following code to 
# first define the hazard function as the p.d.f. divided by the survival function,
weibHaz <- function(x, shape, scale) {
 dweibull(x, shape = shape, scale = scale) / 
  pweibull(x, shape = shape, scale = scale, lower.tail = F)
}

curve(weibHaz(x, shape = 1.5, scale = 1 / 0.03), from = 0, to = 80,
      ylab = "Hazard", xlab = "Time", col = "red", add = TRUE)
curve(weibSurv(x, shape = .75, scale = 1 / 0.75), from = 0, to = 80, ylim = c(.01, .08),
      ylab = "Hazard", xlab = "Time", col = "darkblue", add = TRUE)
curve(weibHaz(x, shape = 1, scale = 1 / 0.03), from = 0, to = 80, 
      ylab = "Hazard", xlab = "Time", col = "black", add = TRUE)

# We may generate random variables from the exponential or Weibull distribution using the 
# functions “rexp” and “rweib”. For example, we may generate 1000 Weibull random variables 
# with shape 1.5 and scale 1/0.03, and compute their mean and median, as follows:
tt.weib <- rweibull(n = 1000, shape = 1.5, scale = 1 / .03)
mean(tt.weib)
median(tt.weib)

# The theoretical mean and median, using equations

# E(T) = (gamma * (1 + 1 / alpha) / lambda)

# and

# t[med] = (log(2) ^ (1 / alpha)) / lambda

# are:
gamma(1 + 1 / 1.5) / .03 # mean

(log(2) ^ (1 / 1.5)) / .03 # median

# The empirical mean and median are close to their theoretical values, as they must be.

# The gamma distribution (!not to be confused with the gamma function!) provides yet another 
# choice for survival modeling. To plot the gamma hazard function for ß = 1.5 and 
# lambda = .03 we can use the following code:
gammaHaz <- function(x, shape, scale) {
 dgamma(x, shape = shape, scale = scale) / 
  pgamma(x, shape = shape, scale = scale, lower.tail = FALSE)
}
curve(gammaHaz(x, shape = 1.5, scale = 1 / .03), from = 0, to = 80,
      ylab = "Hazard", xlab = "Time", col = "red", ylim = c(0, .07))
curve(gammaHaz(x, shape = 1, scale = 1 / .03), from = 0, to = 80,
      ylab = "Hazard", xlab = "Time", col = "black", add = TRUE)
curve(gammaHaz(x, shape = .75, scale = 1 / .03), from = 0, to = 80,
      ylab = "Hazard", xlab = "Time", col = "darkblue", add = TRUE)

# Computing the S(t) from h(t) --------------------------------------------
# If we know the hazard function of a survival random variable, we may derive the survival 
# function using

# µ = E(T) = ∫[∞ - 0] tf(t) dt

# which, when standard integration, and the fact that ƒ(t) = -d / dt * S(t) may be written 
# as

# µ = ∫[∞ - 0] S(t) dt

# For some parametric families, this is simple to do. But if the hazard function is more 
# complicated, we need to use numerical methods to evaluate the integral. We first compute 
# a vector of differences, “tm.diff”, then we find the cumulative hazard functions using the
# “cumsum” function, and finally we use the relationship of the survival function to the 
# cumulative hazard to get the survival outcome of interest.

# End file ----------------------------------------------------------------