# fit_growth.R - Fitting growth curve example

# Reference of generalized logistics:
#   Burger, R et al.: "Comparative analysis of phenomenological growth
#   models applied to epidemic outbreaks", 
#   Mathematical Bioscience and Engineering, 16(5), 4250-4274 (2019).

library(growthrates)
library(doParallel)

grow_logis <- function(time, parms, ...) { # should be equivalence to grow_logistic
  ode_logis <- function(time, init, parms, ...) {
    with(as.list(c(parms, init)), list(mumax * y * (1 - y / K)))
  }
  init <- parms[c("y0")]
  names(init) <- c("y")
  odeparms <- parms[c("mumax", "K")]
  ode(init, time, ode_logis, parms=odeparms, ...)
}

grow_glogis <- function(time, parms, ...) {
  ode_glogis <- function(time, init, parms, ...) {
    with(as.list(c(parms, init)), list(mumax * y^p * (1 - y / K)))
  }
  init <- parms[c("y0")]
  names(init) <- c("y")
  odeparms <- parms[c("mumax", "K", "p")]
  ode(init, time, ode_glogis, parms=odeparms, ...)
}

cAIC <- function(fit) {
  ms <- fit@fit$ms # mean suquared residuals
  n <- nrow(fit@obs)
  k <- length(fit@par)
  # l <- -n / 2 * (log(2 * pi * ms) + 1)
  # AIC <- -2 * l + 2 * k
  # cAIC <- AIC + 2 * k * (k + 1) / (n - k - 1)
  n * (log(2 * pi * ms) + 1) + 2 * k + 2 * k * (k + 1) / (n - k - 1) 
}

fit_growthmodel_multstart <- function(FUN, p, x, y, ...) {
  if (is.data.frame(p))
    p <- as.matrix(p)
  L <- foreach (i=seq_len(nrow(p)), .packages="growthrates") %dopar% {
    try(fit_growthmodel(
      FUN, p[i,], x, y, ...), TRUE)
  }
  attr(L, "convergence") <- sapply(L, function(x) ifelse(
    "try-error" %in% class(x) , NA_real_, x@fit$convergence))
  attr(L, "mse") <- sapply(L, function(x) ifelse(
    "try-error" %in% class(x) , NA_real_, x@fit$ms))
  attr(L, "kappa") <- sapply(L, function(x) ifelse(
    "try-error" %in% class(x) , NA_real_, kappa(x@fit$hessian)))
  L
}

seq.log <- function(from, to, length.out=2) {
  exp(seq(log(from), log(to), length.out=length.out))
}

stop("Expected termination.")

###
y <- c(
  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
  0.000000, 0.100852, 0.100852, 0.201705, 0.201705, 0.201705, 0.201705,
  0.705966, 0.705966, 0.806818, 1.109375, 1.210227, 1.411932, 1.411932,
  1.411932, 1.613636, 2.319602, 2.319602, 3.630682, 4.134943, 4.336648,
  4.437500, 5.244318)
x <- seq_along(y)

r <- diff(y, 2) / diff(x, 2) / y[-c(1, length(y))]
r <- r[!is.na(r) & is.finite(r)]
r <- r[0 < r]

###
lower <- c(y0=1e-5 / 2, mumax=min(r), K=max(y) / 2)
upper <- c(y0=1e5 / 2, mumax=max(r), K=1e5)
p <- c(y0=1e-5, mumax=median(r), K=max(y))
fit <- fit_growthmodel(
  grow_logis, p, x, y, lower=lower, upper=upper, 
  method="Port", control=list(iter.max=300))
summary(fit, cov=FALSE)
print(cAIC(fit))

lower <- c(y0=1e-5 / 2, mumax=min(r), K=max(y) / 2, p=0)
upper <- c(y0=1e5 / 2, mumax=max(r), K=1e5, p=1)
p <- c(y0=1e-5, mumax=median(r), K=max(y), p=1)
fit <- fit_growthmodel(
  grow_glogis, p, x, y, lower=lower, upper=upper,
  method="Port", control=list(iter.max=300, trace=TRUE))
summary(fit, cov=FALSE)
print(cAIC(fit))

###
registerDoParallel(detectCores() / 2)

lower <- c(y0=1e-5 / 2, mumax=min(r), K=max(y) / 2)
upper <- c(y0=1e5 / 2, mumax=max(r), K=1e5)
p <- expand.grid(
  y0=1e-5, mumax=seq.log(min(r) + 1e-6, max(r) - 1e-6, 11), K=seq.log(max(y) / 2 + 1e-6, 1e5 / 2, 11))
L <- fit_growthmodel_multstart(
  grow_logis, p, x, y, lower=lower, upper=upper,
  method="Port", control=list(iter.max=300))
fit <- L[[which.min(attr(L, "mse"))]]
summary(fit, cov=FALSE)
print(cAIC(fit))

lower <- c(y0=1e-5 / 2, mumax=min(r), K=max(y) / 2, p=0)
upper <- c(y0=1e5 / 2, mumax=max(r), K=1e5, p=1)
p <- expand.grid(
  y0=1e-5, mumax=seq.log(min(r) + 1e-6, max(r) - 1e-6, 11), K=seq.log(max(y) / 2 + 1e-6, 1e5 / 2, 11), p=1)
L <- fit_growthmodel_multstart(
  grow_glogis, p, x, y, lower=lower, upper=upper, 
  method="Port", control=list(iter.max=300))
fit <- L[[which.min(attr(L, "mse"))]]
summary(fit, cov=FALSE)
print(cAIC(fit))

###

fit_nCoV_logis <- function(y) {
  x <- seq_along(y)
  
  r <- diff(y, 2) / diff(x, 2) / y[-c(1, length(y))]
  r <- r[!is.na(r) & is.finite(r)]
  r <- r[0 < r]
  if (length(r) <= 0) return(NULL)
  
  lower <- c(y0=1e-5 / 2, mumax=min(r), K=max(y) / 2)
  upper <- c(y0=1e5 / 2, mumax=max(r), K=1e5)
  p <- expand.grid(
    y0=1e-5, mumax=seq.log(min(r) + 1e-6, max(r) - 1e-6, 11), K=max(y))
  
  L <- fit_growthmodel_multstart(
    grow_logis, p, x, y, lower=lower, upper=upper,
    method="Port", control=list(iter.max=300))
  if (all(is.na(attr(L, "mse")))) return(NULL)
  L[[which.min(attr(L, "mse"))]]
}

fit_nCoV_glogis <- function(y) {
  x <- seq_along(y)
  
  r <- diff(y, 2) / diff(x, 2) / y[-c(1, length(y))]
  r <- r[!is.na(r) & is.finite(r)]
  r <- r[0 < r]
  if (length(r) <= 0) return(NULL)
  
  lower <- c(y0=1e-5 / 2, mumax=min(r), K=max(y) / 2, p=0)
  upper <- c(y0=1e5 / 2, mumax=max(r), K=1e5, p=1)
  p <- expand.grid(
    y0=1e-5, mumax=seq.log(min(r) + 1e-6, max(r) - 1e-6, 11), K=max(y), p=1)
    # y0=1e-5, mumax=seq.log(min(r) + 1e-6, max(r) - 1e-6, 11), K=seq.log(lower["K"] + 1e-6, upper["K"] - 1e-6, 11), p=1)

  L <- fit_growthmodel_multstart(
    grow_glogis, p, x, y, lower=lower, upper=upper,
    method="Port", control=list(iter.max=300))
  if (all(is.na(attr(L, "mse")))) return(NULL)
  L[[which.min(attr(L, "mse"))]]
}

stop("Expected termination.")

# t4 is given by cov19.R

L2 <- lapply(t4, fit_nCoV_logis)
AIC2 <- sapply(L2, cAIC)
K2 <- sapply(L2, function(x) x@par)["K",]

L3 <- lapply(t4, fit_nCoV_glogis)
AIC3 <- sapply(L3, cAIC)
K3 <- sapply(L3, function(x) x@par)["K",]

plot(
  AIC2, AIC3, type="n", 
  main="AIC", xlab="logistic", ylab="generalized logistic")
text(AIC2, AIC3, state, col=c(e="red", w="blue", b="purple")[ew])
abline(a=0, b=1, col="gray")

plot(
  K2, K3, type="n", log="xy",
  main="K", xlab="logistic", ylab="generalized logistic")
text(K2, K3, state, col=c(e="red", w="blue", b="purple")[ew])
abline(a=0, b=1, col="gray")
