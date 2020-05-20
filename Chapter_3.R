# =======================================================================================
# Nonparametric Survival Curve Estimation ---------------------------------
# =======================================================================================
library(survival)

# Nonparametric estimation of the Survival Function  ----------------------
# When modeling human or animal survival, it is hard to know what parametric family to 
# choose, and often none of the available families has sufficient flexibility to model the 
# actual shape of the distribution. Thus, in medical and health applications, nonparametric
# methods, which have the flexibility to account for the vagaries of the survival of living
# things, have considerable advantages.

# The most widely used of these is the product-limit estimator, also known as the 
# Kaplan-Meier estimator, which is the product over the failure times of the conditional 
# probabilities of surviving to the next failure time. Formally, it is given by

# S[t] = π[t[i] ≤ t] (1 - q[i]) = π[t[i] ≤ t] (1 - (d[i] / n[i]))

# where n[i] is the number of subjects at risk time t[i], and d[i] is the number of 
# individuals who fail at that time.

# To obtain confidence limits for the product-limit estimator, we first use what is known 
# as the “delta method” to obtain the variance of log(S(t)),

# var(log(S(t))) = ∑[t[i] ≤ t] var(log(1 - q[i])) ≈ ∑[t[i] ≤ t] (d[j] / 
#                                                              n[j] * (n[j] - d[j]))

# To get the variance of S(t) itself, we use the delta method again to obtain

# var(S(t)) ≈ [S(t)] ^ 2 * ∑[t[i] ≤ t] (d[j] / n[j] * (n[j] - d[j]))

# Unfortunately, confidence intervals computed based on this variance may extend above one
# or below zero. While one could truncate them at one and zero, a more satisfying approach
# is to find confidence intervals for the complementary log-log transformation of S(t) as
# follows:

# var(log([-log*S(t)])) ≈ (1 / [S(t)] ^ 2) * ∑[t[i] ≤ t] (d[j] / n[j] * (n[j] - d[j]))

# To obtain estimates of the Kaplan-Meier estimator in R for the following data,
tt <- c(7, 6, 6, 5, 2, 4)
cens <- c(0, 1, 0, 0, 1, 1)
# Note that the “Surv” function produces a special structure for censored survival data:
Surv(tt, cens)

# For the estimation itself we use the “survfit” function. Note that to compute confidence
# intervals based on our preferred method, the complementary log-log transformation, we 
# have to explicitly specify that:
result_KM <- survfit(Surv(tt, cens) ~ 1, conf.type = "log-log")
result_KM
# This prints out the number of “records” (here six patients), the number of patients 
# (n.max and n.start), the number of events (three deaths), the median survival time 
# (6 years), and a 95 % confidence interval for the median. Note that the upper 95 % 
# confidence limit is undefined, indicated by a missing value “NA”. To see the full 
# Kaplan-Meier estimate, and plot it, we use the “summary” and “plot” functions:
summary(result_KM)
plot(result_KM, xlab = "Time to event", ylab = "Survival probability",
     col = c("black", "grey", "grey"))
legend("bottomleft", legend = c("Probability", "95% CI"), 
       lty = c(1, 2), col = c("black", "grey"))

# An alternative estimator of the survival function is known as the Nelson-Altschuler 
# estimator. It is based on the relationship of the survival function to the hazard 
# function. An estimate of the cumulative hazard function is the sum of the estimated 
# hazards up to a time t[i],

# H(t) = ∑[t[i] ≤ t] (d[i] / n[i])

# and the survival function estimate is simply

# S(t) = exp(-H(t))

# We may illustrate this by again referring to the same data as in the afore example. In R, 
# the Nelson-Altschuler estimate may be obtained using the “survfit” function with the
# option “type = "fh", the letters "fh" being taken from the initials of Fleming and 
# Harrington:
result_FH <- survfit(Surv(timeMonths, cens) ~ 1, conf.type = "log-log", type = "fh")
summary(result_FH)

plot(result_FH, xlab = "Time to event", ylab = "Survival probability",
     col = c("black", "grey", "grey"))
legend("bottomleft", legend = c("Probability", "95% CI"), 
       lty = c(1, 2), col = c("black", "grey"))

# =======================================================================================
# Obtaining a Smoothed Hazard and Survival Function -----------------------
# =======================================================================================
# A better way to visualize the hazard function estimate is by using a “kernel” smoother.
# A kernel is a function K.u/, which we center at each failure time. Typically we choose a 
# smooth-shaped kernel, with the amount of smoothing controlled by a parameter b. The 
# estimate of the hazard function is given by:

#  h.hat(t) = 1 / b * ∑[i = 1, D] K * (t - t[i] / b) * d[i] / n[i]

# where t[1] < t[2] < ... < t[D], are distinct ordered failure times.

# While there are many ways to define the kernel function, a common one is the Epanechnikov 
# kernel, K(u) = 3/4 * (1 - µ^2), defined for -1 ≥ µ ≥ 1, and zero elsewhere. In the above 
# formula for the hazard, there is one kernel function placed at each failure time, scaled 
# by the smoothing parameter b. Larger values of b result in wider kernel functions, and 
# hence more smoothing. One problem with this simple approach to hazard estimation is that a 
# kernel may put mass at negative times.

# In the R package, there is a library “muhaz” for estimating and plotting nonparametric 
# hazard functions.
library(muhaz)
t_vec <- c(7, 6, 6, 5, 2, 4)
cens_vec <- c(0, 1, 0, 0, 1, 1)

result_simple <- muhaz(t_vec, cens_vec, max.time = 8, bw.grid = 2.25, bw.method = "global",
                       b.cor = "none")
plot(result_simple)
# The first two arguments are the failure times and censoring indicators, respectively; the 
# maximum time is set at 8; the smoothing parameter b is specified by “bw.grid=2.25”; the 
# “global” option means that a constant smoothing parameter is use for all times; and the 
# “b.cor” option is set to “none” indicating that no boundary correction is to be done.

# =======================================================================================
# Left Truncation ---------------------------------------------------------
# =======================================================================================
# The times between diagnosis and entry into the trial are known as the “backward recurrence
# times.” For patients with times from diagnosis to death (or censoring) is known as 
# “left truncation”. Had a patient died during one of these intervals that patient would not
# have been observed. To obtain an unbiased estimate of the survival distribution, we need 
# to condition on the survival time being greater than the left truncation time. To do this, 
# we construct the Kaplan-Meier estimator as we did earlier, but now a patient only enters 
# the risk set at the left truncation time. Thus, unlike before, the size of the risk set 
# can increase as well as decrease.
tt <- c(7, 6, 6, 5, 2, 4)
status <- c(0, 1, 0, 0, 1, 1)
backTime <- c(-2, -5, -3, -3, -2, -5)
tm_Enter <- -backTime
tm_Exit <- tt - backTime

result_leftTrunc <- survfit(Surv(tm_Enter, tm_Exit, status, type = "counting")
                            ~ 1, conf.type = "none")
summary(result_leftTrunc)
# We have used the terms “tm.enter” and “tm.exit” for the left truncation and survival times, 
# respectively. The reason is derived from the counting process theory, where a subject 
# “enters” the observation period at a particular time and then “exits” it at the time of 
# death or censoring; events that may occur outside of this observation period are not 
# visible to us.

# A serious problem arises with left-truncated data if the risk set becomes empty at an 
# early survival time.


# End file ----------------------------------------------------------------