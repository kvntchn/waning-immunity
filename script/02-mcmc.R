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
	current_frame <- sys.frame(sys.nframe())
	# Solve ODE
	assign("resistant_strain_established", F, envir = current_frame)
	assign("n_to_r", 0, envir = current_frame)
	assign("n_to_wt", 0, envir = current_frame)
	assign("R0", 0, envir = current_frame)
	func <- function(time, state, parameters) {
		compartmental_model(
			time, state, parameters, parent_frame = current_frame
		)}
	withTimeout({
		traj <- data.frame(quietly(ode)(
			y = initial_state,
			times = times,
			parms = theta_proposed,
			func = func,
			method = "impAdams_d"
		)$result)[1:nrow(data),]
	}, timeout = 45 * 2, onTimeout = "warning")
	rm(list = c('n_to_r', 'n_to_wt', 'R0'), envir = current_frame)
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
		lcases <- log(data$cases + 1)
		lincidence <- log(incidence + 1)
		loglik <- dt((lcases - lincidence) / sqrt(
			exp( var(data$cases) - 1 ) * exp( 2 * mean(data$cases) ) + var(data$cases)),
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
		init = parameters,
		constant_which,
		hyperparam = prior_hyperparameters,
		num_iter = 500,
		quiet = F,
		progress = F, acceptance_progress = F,
		C_0,
		batch_size = 50) {
	# Evaluate the function `posterior` at `init`
	post_current <- posterior(theta_proposed = init, hyperparam)
	# Initialize variables to store the current value of theta, the
	# vector of sample values, and the number of accepted proposals
	theta_proposed <- theta_current <- init
	d <- length(init[-constant_which])
	samples <- matrix(NA, ncol = length(init), nrow = num_iter * d)
	colnames(samples) <- names(init)
	accepted <- matrix(F, nrow = num_iter, ncol = d)
	acceptance_rate <- matrix(NA, nrow = num_iter, ncol = d)
	# Run the MCMC algorithm for `num_iter` iterations.
	if (progress)        {pb <- txtProgressBar(min = 0, max = num_iter, style = 3)}
	for (i in 1:num_iter) {
    current_accepted <- rep(0, d)
		# Draw a new theta from a proposal distribution and
		# assign this to a variable called theta_proposed.
		# And update number of accepted proposals.
		for (j in 1:d) {
			theta_proposed[-constant_which][j] <- rnorm(
				1, theta_current[-constant_which][j], exp( C_0[j] * 2 )
			)
			names(theta_proposed) <- names(init)
			# Evaluate the (log) posterior function
			post_proposed <- posterior(theta_proposed, hyperparam)
			log_accept <- (post_proposed - post_current)
			# Compute the Metropolis-Hastings (log) acceptance prob
			if (is.finite(log_accept)) {
				if (log_accept > 0) {
					theta_current <- theta_proposed
					post_current  <- post_proposed
				} else {
					# Draw a random number uniformly-distributed between 0 and 1
					u <- runif(1)
					# Use the random number and the acceptance probability to
					# determine if `theta_proposed` will be accepted.
					if (u < exp(log_accept)) {
						# If accepted, update
						theta_current <- theta_proposed
						post_current  <- post_proposed
						# And update number of accepted proposals.
						current_accepted[j] <- 1
					}}}

			# Add the current theta to the vector of samples.
			output_index <- seq(1, num_iter * d, d)[i] + j - 1
			samples[output_index,] <- theta_current

		}

    accepted[i,] <- as.logical(current_accepted)
		if (i > 1) {
		acceptance_rate[i,] <- colMeans(accepted[1:i,])
		} else {acceptance_rate[i,] <- accepted[1:i,]}
		# Update?
		if (i %% batch_size == 0) {
			which_increase <- which(colMeans(accepted[(i - (batch_size - 1)):i,]) >= 0.455)
			which_decrease <- which(colMeans(accepted[(i - (batch_size - 1)):i,]) <= 0.435)
			C_0[which_increase]  <- C_0[which_increase] + 1 / sqrt(i / batch_size)
			C_0[which_decrease]  <- C_0[which_decrease] - 1 / sqrt(i / batch_size)
		}

		if (!quiet & !progress & i > 1) {
			if (i %% (50/d) == 0 | mean(current_accepted) > 0) {
				# Print the current state of chain and acceptance rate.
				cat("Iteration:", sprintf('% 5d', i),
						"\n\tCurrent:", theta_current[-constant_which],
						"\n\tAcceptance count:",
						paste0(round(sum(rowMeans(accepted[1:i,])))),
						"\n\t Acceptance rate:",
						paste0(round(acceptance_rate[i,] * 100, 2),
									 "% "),
						"\n\n"
				)
			}}
		if (progress) {setTxtProgressBar(pb, i)}
		if (acceptance_progress & i > 1) {
			plot(1:i, cumsum( rowMeans(accepted[1:i,]) ) / (1:i),
					 ylab = "Acceptance rate",
					 ylim = c(0, 0.6), xlim = c(0, i + 10))
			text(i + 8, 0.55,
					 labels = paste0(round(mean(acceptance_rate[i,]) * 100, 2),
					 								"%\n"),
					 pos = 1)
		}
		# i <- i + 1
	}
	return(list(accepted, acceptance_rate, samples, C_0))
}
