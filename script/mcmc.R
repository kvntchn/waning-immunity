# mcmc.R
# February 21, 2023
# Fitting model via MCMC

library(R2jags)
library(coda)

source('script/compartments.R')
source('script/initial_values.R')
san_francisco.dat <- fread("data/san_francisco.csv")
san_francisco.dat <- san_francisco.dat[date >= as.Date('2021-06-15')]
san_francisco.dat[,time := 1:.N]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log prior distributions on parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prior_hyperparameters <- data.frame(
	alpha = c(
		R0_high = 2.2 * 10,
		R0_low = 0.7  * 10,
		disease_duration = 14      * 10,
		recovery_period = 30 * 12  * 10,
		immune_period = 30 * 12    * 10,
		tighten_factor1 = 300  * 1,
		tighten_factor2 = 65  * 1,
		tighten_factor3 = 120  * 1,
		tighten_factor4 = 220  * 1,
		loosen_factor1 = 15  * 1,
		loosen_factor2 = 12  * 1,
		loosen_factor3 = 2  * 1,
		loosen_factor4 = 12  * 1,
	NULL),
	beta = c(
		R0_high = 10,
		R0_low = 10,
		disease_duration = 10,
		recovery_period = 10,
		immune_period   = 10,
		tighten_factor1 = 1,
		tighten_factor2 = 1,
		tighten_factor3 = 1,
		tighten_factor4 = 1,
		loosen_factor1  = 1,
		loosen_factor2 = 1,
		loosen_factor3 = 1,
		loosen_factor4 = 1,
	NULL)
)
prior_hyperparameters$param <- rownames(prior_hyperparameters)

log_prior_theta <- function(theta = parameters, hyperparam = prior_hyperparameters) {
	prior_prob <- sapply(hyperparam$param, function(x) {
		dgamma(theta[[x]],
					 with(hyperparam, alpha[param == x]),
					 with(hyperparam, beta[param == x]),
					 log = T)
	})
	return(sum(unlist(prior_prob)))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log likelihood for trajectory
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_lik_traj <- function(times, data, theta_proposed, initial_state) {
	# Solve ODE
	R0 <<- 0
	traj <- data.frame(quietly(ode)(
		y = initial_state,
		times = times,
		parms = theta_proposed,
		func = compartmental_model,
		method = "lsode"
	)$result)
	rm(list = c('n_to_r', 'n_to_wt', 'R0'), envir = .GlobalEnv)
	# Compute incidence
	beta <- get_beta(traj, theta_proposed)
	incidence <- beta * with(
		traj,
		S * (I_wt + I_r + I_rV) + V * (I_r + I_rV))
	# # Compute vaccinated
	vaccinated <- traj$V * theta_proposed[["vaccination_rate"]]
	return(sum(dpois(data$cases, incidence[1:nrow(data)], log = T)))
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
		log_post(times = seq(from = 1, to = 603, by = theta_proposed[['dt']]),
						 data = san_francisco.dat[,.(cases, deaths, fully_vaccinated)],
						 parameters = theta_proposed,
						 initial_state = get('initial_state', .GlobalEnv)
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
		theta_proposed <- rnorm(
			n = length(theta_current),
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
						"\tCurrent:", theta_current[-which(proposal_sd == 0)],
						"\nAccepted:", c("no", "yes")[current_accepted + 1],
						"\tAcceptance count:",
						paste0(accepted),
						paste0("(", round(accepted / i * 100, 2), "%)\n")
						)
			}}
		if (progress) {setTxtProgressBar(pb, i)}
	}
	return(samples)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN Metropolis-Hastings MCMC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(223)
mcmc_trace <- mh_mcmc(
	posterior = log_post_wrapper,
	init = parameters,
	# init = c(parameters[!names(parameters) %in% c("R0_high", "R0_low")],
					 # parameters[c("R0_high", "R0_low")]),
	proposal_sd = c(rep(0, 6), sqrt(with(prior_hyperparameters, alpha/beta^2)) / 6e4),
	num_iter = 5e5,
	quiet = F,
	progress = F)
save(mcmc_trace, file = 'output/mcmc_trace.rdata')

trace <- mcmc(mcmc_trace[,-(1:8)])

plot(mcmc(mcmc_trace[,7:10]))
