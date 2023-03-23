# mcmc.R
# February 21, 2023
# Fitting model via MCMC

library(R2jags)
library(R.utils)

source('script/01-compartments.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log prior distributions on parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_prior_theta <- function(theta_proposed, hyperparam) {
	if (any(hyperparam$alpha == 0 | hyperparam$beta == 0)) {
		return(-Inf)} else {
			beta.which <- grep("R0_low", hyperparam$param)
			gamma_prior <- sapply(hyperparam$param[-beta.which], function(x) {
				dgamma(transform_parameters(theta_proposed[[x]]),
							 with(hyperparam, alpha[param == x]),
							 with(hyperparam, beta[param == x]),
							 log = T)
			})
			beta_prior <- sapply(hyperparam$param[beta.which], function(x) {
				dbeta(transform_parameters(theta_proposed[[x]]),
							with(hyperparam, alpha[param == x]),
							with(hyperparam, beta[param == x]),
							log = T)
			})
			return(sum(unlist(c(gamma_prior, beta_prior))))
		}
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log likelihood for trajectory
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_lik_traj <- function(times, data = san_francisco.dat, theta_proposed, initial_state) {
	# Solve ODE
	assign("resistant_strain_established", F, envir = .GlobalEnv)
	assign("n_to_r", 0, envir = .GlobalEnv)
	assign("n_to_wt", 0, envir = .GlobalEnv)
	assign("R0", 0, envir = .GlobalEnv)
	withTimeout({
		traj <- data.frame(quietly(ode)(
			y = initial_state,
			times = times,
			parms = theta_proposed,
			func = compartmental_model,
			method = "impAdams_d"
		)$result)[1:nrow(data),]
	}, timeout = 45 * 2, onTimeout = "warning")
	rm(list = c('n_to_r', 'n_to_wt', 'R0'), envir = .GlobalEnv)
	if (exists('traj')) {
		Sys.sleep(0)
		if (nrow(na.exclude(traj)) < nrow(data)) {return(-Inf)}
		N <- sum(initial_state)
		# Compute incidence
		incidence_wt <- with(traj, beta * S * I_wt)
		incidence_r <- with(traj, beta * S * (I_r + I_rV) + V * (I_r + I_rV))
		incidence <- incidence_r + incidence_wt
		incidence[incidence < 0] <- 0
		# Compute vaccinated
		# plot(times, data$cases)
		# lines(times, incidence, col = 'salmon')
		# plot(incidence, data$cases, xlim = c(1, max(data$cases)))
		# abline(coef = c(0,1))
		# loglik <- sum(dpois(round(sum(data$cases)), sum(incidence), log = T))
		# loglik <- c(
		# 	dpois(round(
		# 		sum(data$cases[times < theta_proposed['holiday_date']])),
		# 		sum(incidence[times < theta_proposed['holiday_date']]),
		# 		log = T),
		# 	dpois(round(
		# 		sum(data$cases[times >= theta_proposed['holiday_date'] &
		# 									 	times < theta_proposed['reopening_date']])),
		# 		sum(incidence[times >= theta_proposed['holiday_date'] &
		# 										times < theta_proposed['reopening_date']]),
		# 		log = T),
		# 	dpois(round(
		# 		sum(data$cases[times >= theta_proposed['reopening_date'] &
		# 									 	times < theta_proposed['summer_date']])),
		# 		sum(incidence[times >= theta_proposed['reopening_date'] &
		# 										times < theta_proposed['summer_date']]),
		# 		log = T),
		# 	dpois(round(
		# 		sum(data$cases[times >= theta_proposed['summer_date']])),
		# 		sum(incidence[times >= theta_proposed['summer_date']]),
		# 		log = T)
		# 	)
		# loglik <- dpois(data$cases, incidence, log = T)
		lcases <- log(data$cases + 1)
		lincidence <- log(incidence + 1)
		loglik <- dt((lcases - lincidence) / (sd(data$cases) / mean(data$cases)),
								 length(lcases) - 1,
								 log = T)
		# cat(min(loglik), "\n")
		return( sum( loglik + 7 ) )
	} else {return(-Inf)}
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log posterior
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_post <- function(times, data, theta_proposed, initial_state, hyperparam) {
	log_prior <- log_prior_theta(theta_proposed, hyperparam)
	if (log_prior > -Inf) {
		log_lik   <- log_lik_traj(times, data, theta_proposed, initial_state)
		return(log_prior + log_lik)
	} else {return(-Inf)}
}

# Wrapper
log_post_wrapper <- function(theta_proposed = init, hyperparam) {
	return(
		log_post(times = seq(from = 1, to = nrow(san_francisco.dat),
												 by = theta_proposed[['dt']]),
						 data = san_francisco.dat[,.(cases, deaths, fully_vaccinated)],
						 theta_proposed = theta_proposed,
						 initial_state = get('initial_state', .GlobalEnv),
						 hyperparam
		))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Metropolis-Hastings MCMC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mh_mcmc <- function(
		posterior = log_post_wrapper,
		init = parameters, constant_which,
		hyperparam = prior_hyperparameters,
		num_iter = 500, quiet = F,
		progress = F, acceptance_progress = F,
		s = 2.38, epsilon = 0.05, beta = 0.05,
		C_0,
		adapt_after) {
	# Evaluate the function `posterior` at `init`
	post_current <- posterior(theta_proposed = init, hyperparam)
	# Initialize variables to store the current value of theta, the
	# vector of sample values, and the number of accepted proposals
	theta_current <- init
	d <- length(init[-constant_which])
	samples <- matrix(NA, ncol = length(init), nrow = num_iter)
	colnames(samples) <- names(init)
	accepted <- rep(F, num_iter)
	acceptance_rate <- rep(NA, num_iter)
	# Run the MCMC algorithm for `num_iter` iterations.
	if (progress)        {pb <- txtProgressBar(min = 0, max = num_iter, style = 3)}
	for (i in 1:num_iter) {
		current_accepted <- 0
		# Draw a new theta from a proposal distribution and
		# assign this to a variable called theta_proposed.
		proposal_mean <- theta_current[-constant_which]
		theta_proposed <- c(
			init[constant_which],
			MASS::mvrnorm(
				1,
				proposal_mean,
				Sigma = C_0))
		if (i > adapt_after) {
			# Draw from proposal distribution
			theta_proposed[-constant_which] <- theta_proposed[-constant_which] * beta +
				(1 - beta) * MASS::mvrnorm(
					1,
					proposal_mean,
					Sigma = Rfast::cova(samples[1:(i - 1), -constant_which]))
		}
		names(theta_proposed) <- names(init)
		# Evaluate the (log) posterior function
		post_proposed <- posterior(theta_proposed, hyperparam)
		# Compute the Metropolis-Hastings (log) acceptance prob
		log_accept <- (post_proposed - post_current)
		if (is.finite(log_accept)) {
			if (log_accept > 0) {
				theta_current <- theta_proposed
				post_current <- post_proposed
				# And update number of accepted proposals.
				current_accepted <- 1
			} else {
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
				}}}
		# Add the current theta to the vector of samples.
		samples[i,] <- theta_current
		accepted[i] <- as.logical(current_accepted)
		acceptance_rate[i] <- mean(accepted[1:i])
		if (!quiet & !progress) {
			if (i %% 50 == 0 | current_accepted == 1) {
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
		if (acceptance_progress) {
			plot(1:i, cumsum(accepted[1:i])/(1:i),
					 ylab = "Acceptance rate",
					 ylim = c(0, 0.6), xlim = c(0, i + 10))
			text(i + 8, 0.55, labels = paste0(round(acceptance_rate[i] * 100, 2), "%\n"), pos = 1)
		}
		# i <- i + 1
	}
	return(cbind(accepted, acceptance_rate, samples))
}
