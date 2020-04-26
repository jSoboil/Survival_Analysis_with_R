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
result_FH <- survfit(Surv(tt, cens) ~ 1, conf.type = "log-log", type = "fh")
summary(result_FH)

plot(result_FH, xlab = "Time to event", ylab = "Survival probability",
     col = c("black", "grey", "grey"))
legend("bottomleft", legend = c("Probability", "95% CI"), 
       lty = c(1, 2), col = c("black", "grey"))

# =======================================================================================
# Finding the Median Survival and a Confidence Interval for the Median --------
# =======================================================================================
















