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
	14, 14, 14,
	2, 2, 2, 2,
	0.25, 0.25, 0.25, 0.25,
	7,
	30 * 12, 30 * 12,
	40, 10, 20, 40,
	10, 10, 2, 5) / 2
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

proposal_sd <- c(
	holiday_date =  1,
	reopening_date =  1,
	summer_date =  1,
	R0_high1 = 0.1,
	R0_high2 = 0.1,
	R0_high3 = 0.1,
	R0_high4 = 0.1,
	R0_low1 = 0.01,
	R0_low2 = 0.01,
	R0_low3 = 0.01,
	R0_low4 = 0.01,
	disease_duration = 0.5,
	recovery_period = 1,
	immune_period = 1,
	tighten_factor1 = 1,
	tighten_factor2 = 1,
	tighten_factor3 = 1,
	tighten_factor4 = 1,
	loosen_factor1 = 0.5,
	loosen_factor2 = 0.5,
	loosen_factor3 = 0.1,
	loosen_factor4 = 0.5,
	NULL
)

# # Metropolis in Gibbs
# mcmc_output <- mh_mcmc(
# 	posterior = log_post_wrapper,
# 	init = parameters,
# 	constant_which = constant_which,
# 	num_iter = 2.5e3,
# 	progress = F,
# 	C_0 = log(proposal_sd / 10),
# 	acceptance_progress = T,
# 	batch_size = 100
# 	)
# # save(mcmc_output, file = 'output/mcmc_output_metro-in-gibbs.rdata')

# # Adaptive MH
# mcmc_output <- mh_mcmc(
# 	posterior = log_post_wrapper,
# 	init = parameters,
# 	constant_which = constant_which,
# 	num_iter = 2e5,
# 	progress = T,
# 	C_0 = log(proposal_sd / 2.38),
# 	acceptance_progress = F,
# 	batch_size = 100,
# 	method = "mh"
# 	)
# save(mcmc_output, file = 'output/mcmc_output_mh.rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ODE using parameter estimates from MCMC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load('output/mcmc_output_metro-in-gibbs.rdata')
# load('output/mcmc_output_mh.rdata')
mcmc_trace <- mcmc_output[[3]]
mcmc_trace <- mcmc_trace[
	apply(mcmc_trace[,grep('R0', colnames(mcmc_trace), value = T)], 1, function(x) {all(x > 1e-3)}),]
nrow(mcmc_trace)
plot(mcmc(mcmc_trace[seq(1, nrow(mcmc_trace), 20), -constant_which]))
mcmc_trace_burned <- mcmc_trace[-c(1:(nrow(mcmc_trace)/5)),]
plot(mcmc(mcmc_trace_burned[seq(1, nrow(mcmc_trace_burned), 20), -constant_which]))

mcmc_trace <- as.data.table(mcmc_trace[-(1:(nrow(mcmc_trace)/5)),])
mcmc_parameters_median <- apply(mcmc_trace, 2, quantile, 0.5)
mcmc_parameters_mode <- apply(mcmc_trace, 2, function(x) {x[which.max(density(x)$y)]})

# Initial parameters vs mcmc fit
data.frame(
	prior_mean = parameters,
	posterior_median = mcmc_parameters_median,
	posterior_mode = mcmc_parameters_mode
) -> parameters_estimates
parameters_estimates
# save(parameters_estimates, file = "output/parameters_estimates.rdata")