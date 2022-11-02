# Multivariate state space larval index for capelin 
# load and run capelin_larvae_index.R first

# Load library
library(MARSS)

# 

cap_dat <- cap_reduced %>%  
  group_by(year, station) %>% 
  dplyr::summarise(y = mean(nm2,na.rm=T),
                   d = mean(days_diff,na.rm=T),
                   doy = median(doy, na.rm = T),
                   mean_hatch = mean(mean_hatch, na.rm = T)) %>% droplevels() %>% 
  mutate(station = as.numeric(as.character(station)))

# check autocorrelation
acf(cap_dat$y)
acf(cap_dat$d)
stations <- unique(cap_dat$station)

df <- data.frame(station = stations,
                y = NA)
cap_dat <- full_join(cap_dat,df)

# Use P. Oulette's old strata definitions to at least get an error around estimates
cap_dat %<>% mutate(stratum = ifelse(station %in%
                                               c("8.7", "7.7", "5.7","4.9","3.9","3.8","2.6","1.5","2.5","1.4","2.4","1.3","3.5","3.4",'3.3',"3.2",'3.1',"4.2","4.1","2.3","2.2",'2.1',"1.2","1.1"), "1", 
                                             ifelse(station %in%
                                                      c("10.1","9.4","9.3",'8.2','7.2','7.1',"6.3","6.2","6.1","5.2","5.1",'4.3',"4.4","4.5",'4.6',"3.6","3.7","4.8","5.6",'6.7',"7.5","7.6","6.6"), '2', "3"))) 

cap_dat %<>% group_by(year, stratum) %>% 
  dplyr::summarise(lny = log(mean(y, na.rm = T) + 0.003437986),
                   doy = median(doy, na.rm = T),
                   mean_hatch = mean(mean_hatch, na.rm = T))
                  
cap_dat2 <- cap_dat %>% group_by(year) %>% 
  dplyr::mutate(y = log(mean(y, na.rm = T) + 0.003437986),
                   d = mean(d, na.rm = T)) 

cap_dat %<>% pivot_wider(names_from = stratum, values_from = y)
cap_dat %<>% dplyr::filter(!is.na(year)) 
cap_dat2 %<>% dplyr::filter(!is.na(year)) 
# cap_dat <- cap_dat[rowSums(is.na(cap_dat)) != ncol(cap_dat), ]

# cap_dat <- cap_dat[rowSums(cap_dat[, -(2:4)]) > 0, ]

cap_dat <- as.matrix(cap_dat)
dat <- cap_dat
years = dat[,"year"]
dat = dat[, !(colnames(dat) %in% c("year"))]
dat = t(dat)  #transpose to have years across columns
colnames(dat) = years
n = nrow(dat) - 1

# https://nwfsc-timeseries.github.io/atsa-labs/sec-mss-a-single-well-mixed-population.html Chapter 7.3
# ch 7.3 a single well-mixed population
# assume all stations observations are one pop

# Our model list for a single well-mixed population is:
mod.list.0 <- list(B = matrix(1), U = matrix("u"), Q = matrix("q"), 
                   Z = matrix(1, 3, 1), A = "scaling", R = "diagonal and unequal", 
                   x0 = matrix("mu"), tinitx = 0)
fit.0 <- MARSS(dat, model = mod.list.0)
# Success! abstol and log-log tests passed at 20 iterations.
# Alert: conv.test.slope.tol is 0.5.
# Test with smaller values (<0.1) to ensure convergence.
# 
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Estimation converged in 20 iterations. 
# Log-likelihood: -147.8858 
# AIC: 311.7716   AICc: 313.2716   
# 
# Estimate
# A.2       0.4499
# A.3       1.4191
# R.(1,1)   0.8177
# R.(2,2)   0.3820
# R.(3,3)   0.5502
# U.u       0.0294
# Q.q       0.6490
# x0.mu     1.4445
# Initial states (x0) defined at t=0
# 
# Standard errors have not been calculated. 
# Use MARSSparamCIs to compute CIs and bias estimates.
tidy(fit.0)
glance(fit.0)
plot(fit.0)



# ch 7.4 subpopulations with temporally uncorrelated errors
mod.list.1 <- list(B = "identity", U = "equal", Q = "diagonal and equal", 
                   Z = "identity", A = "scaling", R = "diagonal and unequal", 
                   x0 = "unequal", tinitx = 0)
fit.1 <- MARSS::MARSS(dat, model = mod.list.1)
#Success! algorithm run for 15 iterations. abstol and log-log tests passed.
# Alert: conv.test.slope.tol is 0.5.
# Test with smaller values (<0.1) to ensure convergence.
# 
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Algorithm ran 15 (=minit) iterations and convergence was reached. 
# Log-likelihood: -169.1923 
# AIC: 354.3847   AICc: 355.8847   
# 
# Estimate
# R.(1,1)   0.6117
# R.(2,2)   0.5971
# R.(3,3)   0.9404
# U.1       0.0327
# Q.diag    0.4121
# x0.X.1    0.7944
# x0.X.2    1.3577
# x0.X.3    2.7175
# Initial states (x0) defined at t=0
# 
# Standard errors have not been calculated. 
# Use MARSSparamCIs to compute CIs and bias estimates.
tidy(fit.1)
glance(fit.1)
plot(fit.1)

# ch 7.5 subpopulations with temporally correlated errors
mod.list.2 <- mod.list.1
mod.list.2$Q <- "equalvarcov"
fit.2 <- MARSS::MARSS(dat, model = mod.list.2)
# Success! abstol and log-log tests passed at 348 iterations.
# Alert: conv.test.slope.tol is 0.5.
# Test with smaller values (<0.1) to ensure convergence.
# 
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Estimation converged in 348 iterations. 
# Log-likelihood: -148.1408 
# AIC: 314.2816   AICc: 316.1763   
# 
# Estimate
# R.(1,1)     0.8057
# R.(2,2)     0.3877
# R.(3,3)     0.5488
# U.1         0.0297
# Q.diag      0.6474
# Q.offdiag   0.6468
# x0.X.1      1.3914
# x0.X.2      1.8946
# x0.X.3      2.8949
# Initial states (x0) defined at t=0
# 
# Standard errors have not been calculated. 
# Use MARSSparamCIs to compute CIs and bias estimates.
tidy(fit.2)
glance(fit.2)
plot(fit.2)


par(mfrow = c(2, 2))
for (i in 1:3) {
  plot(years, fit.2$states[i, ], ylab = "log subpopulation estimate", 
       xlab = "", type = "l")
  lines(years, fit.2$states[i, ] - 1.96 * fit.2$states.se[i, 
  ], type = "l", lwd = 1, lty = 2, col = "red")
  lines(years, fit.2$states[i, ] + 1.96 * fit.2$states.se[i, 
  ], type = "l", lwd = 1, lty = 2, col = "red")
  title(rownames(dat)[i])
}


# no strata
mod.list.0 <- list(B = matrix(1), U = matrix("u"), Q = matrix("q"), 
                   Z = matrix(1), A = "scaling", R = "diagonal and unequal", 
                   x0 = matrix("mu"), tinitx = 0)
fit.0 <- MARSS(dat, model = mod.list.0)
# Success! algorithm run for 15 iterations. abstol and log-log tests passed.
# Alert: conv.test.slope.tol is 0.5.
# Test with smaller values (<0.1) to ensure convergence.
# 
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Algorithm ran 15 (=minit) iterations and convergence was reached. 
# Log-likelihood: -53.33989 
# AIC: 114.6798   AICc: 116.0131   
# 
# Estimate
# R.R      0.528
# U.u      0.026
# Q.q      0.414
# x0.mu    2.093
# Initial states (x0) defined at t=0
# 
# Standard errors have not been calculated. 
# Use MARSSparamCIs to compute CIs and bias estimates.
tidy(fit.0)
glance(fit.0)
plot(fit.0)


# Chapter 8, 1 pop with covariate: days diff from spawning 
cap_dat2 <- as.matrix(cap_dat2)
fulldat <- cap_dat2
years = fulldat[, "year"]
fulldat <- fulldat[, !(colnames(fulldat) %in% c("year"))]
dat = fulldat[,1]
dat = t(dat)
colnames(dat) = years
covariates = fulldat[,2]
covariates = t(covariates)
colnames(covariates) = years

# z-score the response variables
the.mean = apply(dat, 1, mean, na.rm = TRUE)
the.sigma = sqrt(apply(dat, 1, var, na.rm = TRUE))
dat = (dat - the.mean) * (1/the.sigma)
# We z-score the covariates to standardize and remove the mean.
the.mean = apply(covariates, 1, mean, na.rm = TRUE)
the.sigma = sqrt(apply(covariates, 1, var, na.rm = TRUE))
covariates = (covariates - the.mean) * (1/the.sigma)

# 8.3 Observation-error only model
# An observation-error only model is a multivariate regression, and we will start here so you see the relationship of MARSS model to more familiar linear regression models.
Q <- U <- x0 <- "zero"
B <- Z <- "identity"
d <- covariates
A <- "zero"
D <- "unconstrained"
y <- dat  # to show relationship between dat & the equation
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, D = D, 
                   d = d, x0 = x0)
kem <- MARSS(y, model = model.list)
# Success! algorithm run for 15 iterations. abstol and log-log tests passed.
# Alert: conv.test.slope.tol is 0.5.
# Test with smaller values (<0.1) to ensure convergence.
# 
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Algorithm ran 15 (=minit) iterations and convergence was reached. 
# Log-likelihood: -46.47204 
# AIC: 96.94409   AICc: 97.30772   
# 
# Estimate
# R.R    0.774
# D.D    0.451
# Initial states (x0) defined at t=0
# 
# Standard errors have not been calculated. 
# Use MARSSparamCIs to compute CIs and bias estimates.
plot(kem)


# 8.4 Process-error only model
# autoregressive process observed without error, and incorporate the covariates into the process model.
R <- A <- U <- "zero"
B <- Z <- "identity"
Q <- "equalvarcov"
C <- "unconstrained"
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, 
                   C = C, c = covariates)
kem <- MARSS(dat, model = model.list)
# Success! algorithm run for 15 iterations. abstol and log-log tests passed.
# Alert: conv.test.slope.tol is 0.5.
# Test with smaller values (<0.1) to ensure convergence.
# 
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Algorithm ran 15 (=minit) iterations and convergence was reached. 
# Log-likelihood: -48.50832 
# AIC: 103.0166   AICc: 103.7666   
# 
# Estimate
# Q.Q      0.867
# x0.x0   -0.979
# C.C      0.139
# Initial states (x0) defined at t=0
# 
# Standard errors have not been calculated. 
# Use MARSSparamCIs to compute CIs and bias estimates.

# better autoregressive model—a mean-reverting model
model.list$B <- "diagonal and unequal"
kem <- MARSS(dat, model = model.list)
# Success! algorithm run for 15 iterations. abstol and log-log tests passed.
# Alert: conv.test.slope.tol is 0.5.
# Test with smaller values (<0.1) to ensure convergence.
# 
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Algorithm ran 15 (=minit) iterations and convergence was reached. 
# Log-likelihood: -42.76039 
# AIC: 93.52077   AICc: 94.8111   
# 
# Estimate
# B.B      0.401
# Q.Q      0.630
# x0.x0   -2.741
# C.C      0.338
# Initial states (x0) defined at t=0
# 
# Standard errors have not been calculated. 
# Use MARSSparamCIs to compute CIs and bias estimates.

# 8.5 Both process- and observation-error
# Here is an example where we have both process and observation error but the covariates only affect the process:
D <- d <- A <- U <- "zero"
Z <- "identity"
B <- "diagonal and unequal"
Q <- "equalvarcov"
C <- "unconstrained"
c <- covariates
# R <- diag(0.16, 2)
R <- A <- U <- "zero"
x0 <- "unequal"
tinitx <- 0
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, 
                   D = D, d = d, C = C, c = c, x0 = x0, tinitx = tinitx)
kem <- MARSS(dat, model = model.list)


# Success! algorithm run for 15 iterations. abstol and log-log tests passed.
# Alert: conv.test.slope.tol is 0.5.
# Test with smaller values (<0.1) to ensure convergence.
# 
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Algorithm ran 15 (=minit) iterations and convergence was reached. 
# Log-likelihood: -42.76039 
# AIC: 93.52077   AICc: 94.8111   
# 
# Estimate
# B.B      0.401
# Q.Q      0.630
# x0.x0   -2.741
# C.C      0.338
# Initial states (x0) defined at t=0
# 
# Standard errors have not been calculated. 
# Use MARSSparamCIs to compute CIs and bias estimates.
tidy(kem)
glance(kem)
plot(kem)
(exp(0.338) - 1) * 100 

resids <- residuals(kem)
acf(resids$.resids)
?acf
# both process and observation error but the covariates only affect the observation process:
C <- c <- A <- U <- "zero"
Z <- "identity"
B <- "diagonal and unequal"
Q <- "equalvarcov"
D <- "unconstrained"
d <- covariates
R <- A <- U <- "zero"
# R <- diag(0.16, 2)
x0 <- "unequal"
tinitx <- 0
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, 
                   D = D, d = d, C = C, c = c, x0 = x0, tinitx = tinitx)
kem <- MARSS(dat, model = model.list, method = "BFGS")
# Success! Converged in 220 iterations.
# Function MARSSkfas used for likelihood calculation.
# 
# MARSS fit is
# Estimation method: BFGS 
# Estimation converged in 220 iterations. 
# Log-likelihood: -43.76338 
# AIC: 95.52677   AICc: 96.81709   
# 
# Estimate
# B.B      0.439
# Q.Q      0.666
# x0.x0   -2.375
# D.D      0.244
# Initial states (x0) defined at t=0
# 
# Standard errors have not been calculated. 
# Use MARSSparamCIs to compute CIs and bias estimates.
tidy(kem)
glance(kem)

# a way to include multiple strata with this?? east west?
df<-as.data.frame(t(kem$states))
df2 <- as.data.frame(t(kem$states.se))
