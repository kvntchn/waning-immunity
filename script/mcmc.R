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
log_prior_theta <- function(theta = parameters) {
	R0_high <- dunif(
		theta[["R0_high"]], 1, 3, log = T)
	R0_low  <- dunif(
		theta[["R0_low"]], 0.375, 1.5, log = T)
	R0_omicron  <- dunif(
		theta[["R0_omicron"]], 1.5, 6, log = T)
	disease_duration <- dunif(
		theta[["disease_duration"]], 14, 28, log = T)
	disease_duration_omicron <- dunif(
		theta[["disease_duration_omicron"]], 7, 21, log = T)
	death_rate <- dunif(
		theta[['death_rate']], 1e-4, 1e-2, log = T)
	vaccination_rate <- dunif(
		theta[['vaccination_rate']], 1e-4, 1e-2, log = T)
	recovery_period <- dunif(
		theta[['recovery_period']], 30 * 6, 30 * 12 * 1.5, log = T)
	recovery_period_omicron <- dunif(
		theta[['recovery_period_omicron']], 30 * 6, 30 * 12, log = T)
	emergence_probability <- dunif(
		theta[['emergence_probability']], 1e-5, 1e-3, log = T)
	omicron_date <- dunif(theta[['omicron_date']], 360, 420, log = T)
	return(R0_high + R0_low + R0_omicron + disease_duration +
				 	disease_duration_omicron + death_rate +
				 	vaccination_rate + recovery_period +
				 	recovery_period_omicron + emergence_probability +
				 	omicron_date
	)}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log likelihood for single observation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_lik_I <- function(data, model) {
	# Incidence is observed through a Poisson process.
	return(dpois(
		x = data,
		lambda = model,
		log = TRUE
	))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log likelihood for entire trajectory
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_lik_traj <- function(times, data, theta_proposed, initial_state) {
	traj <- data.frame(ode(
		y = initial_state,
		times = times,
		parms = theta_proposed,
		func = compartmental_model,
		method = "lsode"
	))
	model <- traj$I_wt + traj$I_r + traj$I_rV
	return(sum(log_lik_I(data, model)))
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
		if (abs(log_lik) < .Machine$double.eps) {
			return(-Inf)} else {
				return(log_prior + log_lik)}
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
mh_mcmc <- function(posterior, init, proposal_sd, num_iter) {
	# Evaluate the function `posterior` at `init`
	post_current <- posterior(init)
	# Initialize variables to store the current value of theta, the
	# vector of sample values, and the number of accepted proposals
	theta_current <- init
	samples <- matrix(NA, ncol = length(init), nrow = num_iter + 1)
	samples[1,] <- init
	accepted <- 0
	# Run the MCMC algorithm for `num_iter` iterations.
	for (i in 1:num_iter) {
		# Draw a new theta from a Gaussian proposal distribution and
		# assign this to a variable called theta_proposed.
		theta_proposed <- rnorm(n = length(theta_current),
														mean = theta_current,
														sd = proposal_sd)
		names(theta_proposed) <- names(theta_current)
		# Evaluate the (log) posterior function
		post_proposed <- posterior(theta_proposed)
		# Compute the Metropolis-Hastings (log) acceptance rate
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
				accepted <- accepted + 1
			}}
		# Add the current theta to the vector of samples.
		samples[i,] <- theta_current
		# Print the current state of chain and acceptance rate.
		cat("Iteration:", i, "Chain:", theta_current,
				"acceptance rate:", accepted / i, "\n")
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
	p_non_vax = 0.2,
	saturation = 0.01,
	emergence_probability = 5e-4,
	omicron_date = 390)

mcmc_trace <- mh_mcmc(
	posterior = log_post_wrapper,
	init = parameters,
	proposal_sd = c(
		rep(0, 4),
		rep(0.001, length(parameters) - 4)
	),
	num_iter = 20)


compartmental_model <- function() {

	traj <- data.frame(ode(
		y = initial_state,
		times = times,
		parms = theta_proposed,
		func = compartmental_model,
		method = "lsode"
	))
	model <- traj$I_wt + traj$I_r + traj$I_rV

  for (k in 1:K) {
  	I[k] ~ dpois(lambda[k])
  	lambda[k]
  }

    logit(pi[i]) <- b[1] + b[2]*tobgp2[i] + b[3]*tobgp3[i] + b[4]*tobgp4[i] + b[5]*age[i];
    case[i] ~ dbin(pi[i],1);
  }

  # PRIORS ON BETAS
  for (j in 1:Nx) {b[j] ~ dnorm(mu[j],tau[j]);} # independent normal prior

  # Calculate ORs:
  for (l in 1:Nx) {
      OR[l] <- exp(b[l]);
  }
}
