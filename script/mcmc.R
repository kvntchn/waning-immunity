# mcmc.R
# February 21, 2023
# Fitting model via MCMC

library(R2jags)
library(coda)

source('script/compartments.R')
san_francisco.dat <- fread("data/san_francisco.csv")
san_francisco.dat <- san_francisco.dat[-(1:75),]
san_francisco.dat[,time := 1:.N]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log prior distributions on parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prior_hyperparameters <- data.frame(
# 	param = c("R0_high", "R0_low", "R0_omicron",
# 						"disease_duration", "disease_duration_omicron",
# 						"death_rate", "vaccination_rate",
# 						"recovery_period", "recovery_period_omicron",
# 						"emergence_probability",
# 						"omicron_date"),
# 	alpha = c(1.5, 7.5, 30,
# 						170, 150,
# 						20, 20,
# 						30 * 12 / 10, 30 * 10 / 10,
# 						5e-4 * 1e-2,
# 						390 * 3
# 	),
# 	beta = c(10, 10, 10,
# 					 10, 10,
# 					 20000, 20000,
# 					 1 / 10, 1 / 10,
# 					 1 * 1e-2, 3)
# )

prior_hyperparameters <- data.frame(
	param = c("R0_high", "R0_low", "R0_omicron",
						"disease_duration", "disease_duration_omicron",
						"death_rate", "vaccination_rate",
						"recovery_period", "recovery_period_omicron",
						"emergence_probability",
						"omicron_date"),
	alpha = c(1, 1/100, 1,
						3, 3,
						1e-6, 1e-6,
						1, 1,
						1e-6,
						390
	),
	beta = c(100, 1, 100,
					 300, 300,
					 1e-1, 1e-1,
					 30 * 12 * 10, 30 * 10 * 10,
					 5e-2,
					 410)
)

log_prior_theta <- function(theta = parameters, hyperparam = prior_hyperparameters) {
	prior_prob <- sapply(names(theta), function(x) {
		dunif(theta[[x]],
					with(hyperparam, alpha[param == x]),
					with(hyperparam, beta[param == x]),
					log = T)
	})
	return(with(prior_prob,
							R0_high + R0_low + R0_omicron +
								disease_duration + disease_duration_omicron +
								death_rate + vaccination_rate +
								recovery_period + recovery_period_omicron +
								emergence_probability + omicron_date
	)
	)}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log likelihood for trajectory
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_lik_traj <- function(times, data, theta_proposed, initial_state) {
	traj <- data.frame(quietly(ode)(
		y = initial_state,
		times = times,
		parms = theta_proposed,
		func = compartmental_model,
		method = "lsode"
	)$result)
	model <- traj$I_wt + traj$I_r + traj$I_rV
	return(sum(dpois(data, model[1:length(data)], log = T)))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log posterior
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_post <- function(times, data, parameters, initial_state) {
	log_prior <- log_prior_theta(parameters)
	if (!is.finite(log_prior)) {
		return(log_prior)
	} else {
		log_lik <- log_lik_traj(times, data, parameters, initial_state)
		if (!is.finite(log_lik)) {
			return(-Inf)
		} else {
			if (abs(log_lik) < .Machine$double.eps) {
				return(-Inf)
			} else {
				return(log_prior + log_lik)
			}
		}
	}
}

# Wrapper
log_post_wrapper <- function(theta_proposed) {
	return(
		log_post(times = seq(from = 0, to = 930, by = theta_proposed[['dt']]),
						 data = san_francisco.dat$cases,
						 parameters = theta_proposed,
						 initial_state = c(
						 	S = N - 30,
						 	I_wt = 30,
						 	I_r = 0,
						 	I_rV = 0,
						 	R = 0,
						 	RV = 0,
						 	D = 0,
						 	V = 0)
		))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Metropolis-Hastings MCMC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mh_mcmc <- function(posterior, init, proposal_sd, num_iter, quiet = T, progress = F) {
	# Evaluate the function `posterior` at `init`
	post_current <- posterior(init)
	# Initialize variables to store the current value of theta, the
	# vector of sample values, and the number of accepted proposals
	theta_current <- init
	samples <- matrix(NA, ncol = length(init), nrow = num_iter + 1)
	colnames(samples) <- names(init)
	samples[1,] <- init
	accepted <- 0
	# Run the MCMC algorithm for `num_iter` iterations.
	if (progress) {pb <- txtProgressBar(min = 0, max = num_iter, style = 3)}
	for (i in 2:(num_iter + 1)) {
		current_accepted <- 0
		# Draw a new theta from a Gaussian proposal distribution and
		# assign this to a variable called theta_proposed.
		theta_proposed <- rnorm(n = length(theta_current),
														mean = theta_current,
														sd = proposal_sd)
		names(theta_proposed) <- names(init)
		# Evaluate the (log) posterior function
		post_proposed <- posterior(theta_proposed)
		# Compute the Metropolis-Hastings (log) acceptance rate
		log_accept <- post_proposed - post_current
		if (is.finite(log_accept)) {
			# Draw a random number uniformly-distributed between 0 and 1
			u <- runif(1)
			# Use the random number and the acceptance probability to
			# determine if `theta_proposed` will be accepted.
			if (log(u) < log_accept) {
				# If accepted, update
				theta_current <- theta_proposed
				post_current <- post_proposed
				# And update number of accepted proposals.
				accepted <- accepted + 1
				current_accepted <- 1
			}}
		# Add the current theta to the vector of samples.
		samples[i,] <- theta_current
		if (!quiet) {
			if (i %% 500 == 0 | current_accepted == 1) {
				# Print the current state of chain and acceptance rate.
				cat("Iteration:", sprintf('% 5d', i),
						"\tChain:", theta_current[-which(proposal_sd == 0)],
						"\nAccepted:", c("no", "yes")[current_accepted + 1],
						"\tAcceptance rate:", round(accepted / i * 100, 2), "%\n")
			}}
		if (progress) {setTxtProgressBar(pb, i)}
	}
	return(samples)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN Metropolis-Hastings MCMC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N <- 8e5
parameters <- c(
	# Fixed parameters
	F_h = N / 2000,
	R0_previous = -Inf,
	immune_period = Inf,
	dt = 1,
	p_non_vax = 0.2,
	saturation = 0.01,
	# Random parameters
	R0_high = 1.5,
	R0_low = 0.75,
	R0_omicron = 3,
	disease_duration = 17,
	disease_duration_omicron = 15,
	death_rate = 1 / 1000,
	vaccination_rate = 1 / 1000,
	recovery_period = 30 * 12,
	recovery_period_omicron = 30 * 10,
	emergence_probability = 5e-4,
	omicron_date = 390)

init <- parameters

set.seed(230)
mcmc_trace <- mh_mcmc(
	posterior = log_post_wrapper,
	init = init,
	proposal_sd = c(
		rep(0, 6),
		sqrt((prior_hyperparameters$beta - prior_hyperparameters$alpha)^2 / 12) / 100
	),
	num_iter = 5e4,
	quiet = F,
	progress = F)
save(mcmc_trace, file = 'output/mcmc_trace.rdata')

trace <- mcmc(mcmc_trace[,-(1:6)])

plot(mcmc(mcmc_trace[,7:9]))
