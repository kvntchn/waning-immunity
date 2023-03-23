# mcmc.R
# February 21, 2023
# Fitting model via MCMC

library(data.table)
library(R2jags)
san_francisco.dat <- fread("data/san_francisco.csv")
san_francisco.dat <- san_francisco.dat[date >= as.Date('2021-06-15')]
san_francisco.dat[,time := 1:.N]

source("script/02-mcmc.R")
source('script/03-initial_values.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN Metropolis-Hastings MCMC (No emergence)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
transform_parameters <- function(parameters) {
	parameters[grep("R0_high|loosen_factor|tighten_factor", names(parameters))] <-
		parameters[grep("R0_high|loosen_factor|tighten_factor", names(parameters))] - 1
	return(parameters)
}
constant_which <- 1:grep("emergence_date", names(parameters))
d <- length(parameters[-constant_which])
prior_mean <- transform_parameters(parameters[-constant_which])
prior_sd   <- c(
	5, 5, 5,
	0.5, 0.5, 0.5, 0.5,
	0.05, 0.05, 0.05, 0.05,
	7,
	30 * 6, 30 * 6,
	80, 20, 80, 80,
	4, 4, 1, 2)
prior_hyperparameters <- data.frame(
	param = names(prior_mean),
	alpha = prior_mean^2 / prior_sd^2,
	beta  = prior_mean / prior_sd^2
)

R0_low.which <- grep("R0_low", names(prior_mean))
prior_hyperparameters[R0_low.which,'alpha'] <- (
	((1 - prior_mean) / prior_sd^2 - 1 / prior_mean) * prior_mean^2)[R0_low.which]
prior_hyperparameters[R0_low.which,'beta'] <- (
	prior_hyperparameters$alpha * (1 / prior_mean - 1)
)[R0_low.which]

# ggplot() +
# 	geom_function(fun = function(x) {
# 		with(prior_hyperparameters, dgamma(x - 1, alpha[15], beta[15]))
# 	}, aes(col = '1')) +
# 	geom_function(fun = function(x) {
# 		with(prior_hyperparameters, dgamma(x - 1, alpha[16], beta[16]))
# 	}, aes(col = '2')) +
# 	geom_function(fun = function(x) {
# 		with(prior_hyperparameters, dgamma(x - 1, alpha[17], beta[17]))
# 	}, aes(col = '3')) +
# 	geom_function(fun = function(x) {
# 		with(prior_hyperparameters, dgamma(x - 1, alpha[18], beta[18]))
# 	}, aes(col = '4')) +
# 	scale_x_continuous(limits = c(1, 500)) +
# 	theme_bw()

# prior_uniform_hyperparameters <- data.frame(
# 	param = names(prior_mean),
# 	alpha = c(#150, 250, 375,
# 						1.5, 2, 1.5, 1.5,
# 						0, 0, 0, 0,
# 						0, 1, 1,
# 						100, 10, 100, 100,
# 						1, 1, 1, 1,
# 						NULL
# 						),
# 	beta  = c(#180, 320, 420,
# 						5, 15, 5, 5,
# 						1, 1, 1, 1,
# 						30, 30 * 12 * 2, 30 * 12 * 2,
# 						500, 400, 400, 500,
# 						50, 200, 50, 50,
# 						NULL
# 						)
# )

set.seed(222)

mcmc_trace_no_emergence <- mh_mcmc(
	posterior = log_post_wrapper,
	init = parameters,
	constant_which = constant_which,
	num_iter = 2.5e4,
	progress = F,
	C_0 = diag((prior_sd * 0.1)^2 / d, d),
	adapt_after = d * 2,
	epsilon = 0.05,
	beta = 0.05
	)
# save(mcmc_trace_no_emergence, file = 'output/mcmc_trace_no_emergence.rdata')

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ODE using parameter estimates from MCMC
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # load('output/mcmc_trace_no_emergence.rdata')
# plot(mcmc(mcmc_trace_no_emergence[seq(1, nrow(mcmc_trace_no_emergence), 40),
# 																	-(c(constant_which, 1:2 + max(constant_which)))]))
# mcmc_trace_no_emergence_burned <- mcmc_trace_no_emergence[-(1:nrow(mcmc_trace_no_emergence)/2),]
# plot(mcmc(mcmc_trace_no_emergence_burned[seq(1, nrow(mcmc_trace_no_emergence_burned), 5), -(1:13)]))
# mcmc_trace <- as.data.table(mcmc_trace_no_emergence[-(1:3e4),])
# mcmc_parameters <- colMeans(mcmc_trace)
#
# # Initial parameters vs mcmc fit
# data.frame(
# 	prior_mean = parameters,
# 	fitted = mcmc_parameters
# )
#
