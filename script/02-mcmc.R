# mcmc.R
# February 21, 2023
# Fitting model via MCMC

library(R2jags)

source('script/01-compartments.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log prior distributions on parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_prior_theta <- function(theta = parameters, hyperparam) {
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
	assign("R0", 0, envir = .GlobalEnv)
	traj <- data.frame(quietly(ode)(
		y = initial_state,
		times = times,
		parms = theta_proposed,
		func = compartmental_model,
		method = "lsode"
	)$result)
	rm(list = c('n_to_r', 'n_to_wt', 'R0'), envir = .GlobalEnv)
	N <- sum(initial_state)
	# Compute incidence
	beta <- get_beta(traj, theta_proposed)
	incidence <- beta * with(
		traj,
		S * (I_wt + I_r + I_rV) + V * (I_r + I_rV))
	# Compute vaccinated
	# vaccinated <- traj$V * theta_proposed[["vaccination_rate"]]
	# loglik <- dpois(sum(data$cases), sum(incidence[1:nrow(data)]), log = T)
	# loglik <- sum( dpois( data$cases, incidence[1:nrow(data)], log = T) )
	loglik <- dbinom(sum(data$cases), size = N, prob = sum(incidence[1:nrow(data)]) / N, log = T)
	return(loglik)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log posterior
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_post <- function(times, data, parameters, initial_state, hyperparam) {
	log_prior <- log_prior_theta(parameters, hyperparam)
	log_lik   <- log_lik_traj(times, data, parameters, initial_state)
	return(log_prior + log_lik)
}

# Wrapper
log_post_wrapper <- function(theta_proposed, hyperparam) {
	return(
		log_post(times = seq(from = 1, to = nrow(san_francisco.dat),
												 by = theta_proposed[['dt']]),
						 data = san_francisco.dat[,.(cases, deaths, fully_vaccinated)],
						 parameters = theta_proposed,
						 initial_state = get('initial_state', .GlobalEnv),
						 hyperparam
		))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Metropolis-Hastings MCMC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mh_mcmc <- function(
		posterior, init, constant_which,
		hyperparam = prior_hyperparameters,
		num_iter, quiet = F,
		progress = F, acceptance_progress = F) {
	# Evaluate the function `posterior` at `init`
	post_current <- posterior(init, hyperparam)
	# Initialize variables to store the current value of theta, the
	# vector of sample values, and the number of accepted proposals
	theta_current <- init
	samples <- matrix(NA, ncol = length(init), nrow = num_iter)
	colnames(samples) <- names(init)
	accepted <- rep(NA, num_iter)
	acceptance_rate <- rep(NA, num_iter)
	# Run the MCMC algorithm for `num_iter` iterations.
	if (progress)        {pb <- txtProgressBar(min = 0, max = num_iter, style = 3)}
	if (acceptance_progress) {arb <- txtProgressBar(min = 0, max = 1, style = 3)}
	for (i in 1:(num_iter)) {
		current_accepted <- 0
		# Draw a new theta from a proposal distribution and
		# assign this to a variable called theta_proposed.
		d <- length(theta_current[-constant_which])
		proposal_mean <- theta_current[-constant_which]
		theta_proposed <- c(
			init[constant_which],
			rnorm(
				n = d,
				proposal_mean,
				0.1 / sqrt(d))
		)
		# http://www2.stat.duke.edu/~scs/Courses/Stat376/Papers/AdaptiveMC/RobertsRosenthalAdaptExamples2006.pdf
		if (i > 2 * d) {
			theta_proposed[-constant_which] <- theta_proposed[-constant_which]  * 0.05 +
				(1 - 0.05) * MASS::mvrnorm(
					1,
					proposal_mean,
					2.38^2 / (d) * cov(samples[1:(i - 1),-constant_which]))
		}
		names(theta_proposed) <- names(init)
		# Evaluate the (log) posterior function
		post_proposed <- posterior(theta_proposed, hyperparam)
		# Compute the Metropolis-Hastings (log) acceptance prob
		log_accept <- post_proposed - post_current
		if (is.finite(log_accept)) {
			# Draw a random number uniformly-distributed between 0 and 1
			u <- runif(1)
			# Use the random number and the acceptance probability to
			# determine if `theta_proposed` will be accepted.
			if (u < exp(log_accept)) {
				# If accepted, update
				theta_current <- theta_proposed
				post_current <- post_proposed
				# And update number of accepted proposals.
				current_accepted <- 1
			}}
		# Add the current theta to the vector of samples.
		samples[i,] <- theta_current
		accepted[i] <- current_accepted
		acceptance_rate[i] <- mean(accepted[1:i])
		if (!quiet & !progress) {
			if (i %% 500 == 0 | current_accepted == 1) {
				# Print the current state of chain and acceptance rate.
				cat("Iteration:", sprintf('% 5d', i),
						"\tCurrent:", theta_current[-constant_which],
						"\nAccepted:", c("no", "yes")[current_accepted + 1],
						"\tAcceptance count:",
						paste0(sum(accepted[1:i])),
						paste0("(", round(acceptance_rate[i] * 100, 2), "%)\n")
				)
			}}
		if (progress) {setTxtProgressBar(pb, i)}
		if (acceptance_progress) {setTxtProgressBar(arb, acceptance_rate)}
	}
	return(cbind(accepted, acceptance_rate, samples))
}
