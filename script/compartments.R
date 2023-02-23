# state_model.R
# February 9, 2023
# Simple ODE implementation

library(tidyverse)
library(data.table)
library(deSolve)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Derivatives function for closed compartmental model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compartmental_model <- function(time, state = initial_state, parameters) {
	# Parameters:
	F_h     <- parameters['F_h']
	gamma   <- 1 / parameters["disease_duration"]
	delta   <- parameters["death_rate"]
	theta0  <- parameters["vaccination_rate"]
	w       <- 1 / parameters["immune_period"]
	mu      <- 1 / parameters["recovery_period"]
	p       <- parameters['emergence_probability']
	R0_high <- parameters["R0_high"]
	R0_low  <- parameters["R0_low"]
	tighten_factor1   <- parameters["tighten_factor1"]
	tighten_factor2   <- parameters["tighten_factor2"]
	tighten_factor3   <- parameters["tighten_factor3"]
	tighten_factor4   <- parameters["tighten_factor4"]
	loosen_factor1   <- parameters["loosen_factor1"]
	loosen_factor2   <- parameters["loosen_factor2"]
	loosen_factor3   <- parameters["loosen_factor3"]
	loosen_factor4   <- parameters["loosen_factor4"]
	h       <- parameters["p_non_vax"]
	k       <- parameters["saturation"]
	dt      <- parameters["dt"]
	R0     <- get('R0', envir = .GlobalEnv)
	# States:
	S    <- state["S"]
	I_wt <- state["I_wt"]
	I_r  <- state["I_r"]
	I_rV <- state["I_rV"]
	R    <- state["R"]
	RV   <- state["RV"]
	D    <- state["D"]
	V    <- state["V"]
	N <-  S + I_wt + I_r + I_rV + R + RV + D + V
	# Time-varying parameters
	## Vaccination rate
	theta <- (1 - h / (S + R + I_wt + I_r)) * theta0 / (S + R + k)
	## Time-varying force of  infection
	if (time > 0)   {F_h <- N / tighten_factor1; loosen_factor <- loosen_factor1}
	# time 160 corresponds to 11/21/2021 (Thanksgiving)
	if (time >= 160) {F_h <- N / tighten_factor2; loosen_factor <- loosen_factor2}
	# time 305 corresponds to 04/15/2022 (SF reopening)
	if (time >= 305) {F_h <- N / tighten_factor3; loosen_factor <- loosen_factor3}
	# time 321 corresponds to 05/01/2022 (Warmer weather)
	if (time >= 410) {F_h <- N / tighten_factor4; loosen_factor <- loosen_factor4}
	if ((I_wt + I_r + I_rV) > F_h & R0 > R0_low)                   {R0 <<- R0_low}
	if ((I_wt + I_r + I_rV) < (F_h / loosen_factor) & R0 < R0_high)  {R0 <<- R0_high}
	beta  <- get('R0', envir = .GlobalEnv) * gamma / N
	# Resistant strain?
	if ((I_r + I_rV) <= 0.0001 * N) {
		n_to_r  <<- rpois(1, dt * p * max(0, I_wt, na.rm = T))
		n_to_wt <<- rpois(1, dt * p * max(0, I_r, na.rm = T))
	I_r  <- I_r  - n_to_wt + n_to_r
	I_wt <- I_wt + n_to_wt - n_to_r
	}
	# Derivatives:
	dS    <- mu * R + w * V - theta * S - beta * (I_wt + I_r + I_rV) * S
	dI_wt <- - (gamma + delta) * I_wt + beta * S * (I_wt)
	dI_r  <- - (gamma + delta) * I_r  + beta * S * (I_r + I_rV) + p * I_wt
	dI_rV <- - (gamma + delta) * I_rV + beta * V * (I_r + I_rV)
	dR    <- - mu * R - theta * R + gamma * (I_wt + I_r)
	dRV   <- - mu * RV + theta * R + gamma * I_rV
	dD    <- delta * (I_wt + I_r + I_rV)
	dV    <- mu * RV + theta * S - beta * V + (I_r + I_rV) - w * V
	return(list(c(dS, dI_wt, dI_r, dI_rV, dR, dRV, dD, dV)))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute what force of infection should be given state
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_beta <- function(traj, parameters) {
	beta <- rep(NA, nrow(traj))
	tighten_factor <- beta
	tighten_factor[traj[,'time'] > 0]    <- parameters['tighten_factor1']
	tighten_factor[traj[,'time'] >= 160] <- parameters['tighten_factor2']
	tighten_factor[traj[,'time'] >= 305] <- parameters['tighten_factor3']
	tighten_factor[traj[,'time'] >= 410] <- parameters['tighten_factor4']
	loosen_factor <- beta
	loosen_factor[traj[,'time'] > 0]    <- parameters['loosen_factor1']
	loosen_factor[traj[,'time'] >= 160] <- parameters['loosen_factor2']
	loosen_factor[traj[,'time'] >= 305] <- parameters['loosen_factor3']
	loosen_factor[traj[,'time'] >= 410] <- parameters['loosen_factor4']
	R0 <- c(beta, NA)
	R0[1] <- 0
	for (i in 1:length(beta)) {
		I <- (traj[i, 'I_wt'] + traj[i, 'I_r'] + traj[i, 'I_rV'])
		R0[i + 1] <- R0[i]
		if (I > N / tighten_factor[i] & R0[i] > parameters['R0_low']) {
			R0[i + 1] <- parameters['R0_low']}
		if (I < (N / tighten_factor[i] / loosen_factor[i]) & R0[i] < parameters['R0_high'])  {
			R0[i + 1] <- parameters['R0_high']}
	}
	return(R0[-1] / parameters['disease_duration'] / N)
}