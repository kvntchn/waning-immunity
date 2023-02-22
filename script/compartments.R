# state_model.R
# February 9, 2023
# Simple ODE implementation

library(tidyverse)
library(data.table)
library(deSolve)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Derivatives function for closed compartmental model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compartmental_model <- function(time, state, parameters) {
	# Parameters:
	F_h     <- parameters['F_h']
	gamma_delta   <- 1 / parameters["disease_duration"]
	gamma_omicron <- 1 / parameters["disease_duration_omicron"]
	delta   <- parameters["death_rate"]
	theta0  <- parameters["vaccination_rate"]
	w       <- 1 / parameters["immune_period"]
	mu_delta      <- 1 / parameters["recovery_period"]
	mu_omicron      <- 1 / parameters["recovery_period_omicron"]
	p       <- parameters['emergence_probability']
	R0_high_delta <- parameters["R0_high"]
	R0_low  <- parameters["R0_low"]
	R0_omicron  <- parameters["R0_omicron"]
	if (time <= 1) {R0_prev <<- parameters["R0_previous"]}
	h       <- parameters["p_non_vax"]
	k       <- parameters["saturation"]
	omicron_date <- parameters["omicron_date"]
	dt      <- parameters["dt"]
	# States:
	S    <- state["S"]
	I_wt <- state["I_wt"]
	I_r  <- state["I_r"]
	I_rV <- state["I_rV"]
	R    <- state["R"]
	RV   <- state["RV"]
	D    <- state["D"]
	V    <- state["V"]
	N <- S + I_wt + I_r + I_rV + R + RV + D + V
	I <- I_wt + I_r + I_rV
	# Time-varying parameters
	## Duration of disease
	if (time < omicron_date) {
		R0_high <<- R0_high_delta
		gamma_current <<- gamma_delta
	} else {
		R0_high <<- R0_omicron
		gamma_current <<- gamma_omicron
	}
	## Recovery period
	if (time < omicron_date) {
		mu <<- mu_delta
	} else {
		mu <<- mu_omicron
	}
	## Vaccination rate
	theta <- (1 - h / (S + R + I_wt + I_r)) * theta0 / (S + R + k)
	## Time-varying force of  infection
	if (I > F_h       & R0_prev > R0_low)  {R0 <<- R0_low}
	if (I < (F_h / 8) & R0_prev < R0_high) {R0 <<- R0_high}
	## Holiday period
	if (time >= 415 & time <= 425) {R0 <<- R0_high + 3}
	## Back to bar rush
	if (time >= 569 & time <= 579) {R0 <<- R0_high + 1}
	beta  <- R0 * gamma_current / N
	R0_prev <<- R0
	# Resistant strain?
	if ((I_r + I_rV) <= 0.0001 * N) {
		n_to_r  <<- rpois(1, dt * p * max(0, I_wt, na.rm = T))
		n_to_wt <<- rpois(1, dt * p * max(0, I_r, na.rm = T))
	I_r  <- I_r  - n_to_wt + n_to_r
	I_wt <- I_wt + n_to_wt - n_to_r
	}
	# Derivatives:
	dS    <- mu * R + w * V - theta * S - beta * (I_wt + I_r + I_rV) * S
	dI_wt <- - (gamma_current + delta) * I_wt + beta * S * (I_wt)
	dI_r  <- - (gamma_current + delta) * I_r  + beta * S * (I_r + I_rV) + p * I_wt
	dI_rV <- - (gamma_current + delta) * I_rV + beta * V * (I_r + I_rV)
	dR    <- - mu * R - theta * R + gamma_current * (I_wt + I_r)
	dRV   <- - mu * RV + theta * R + gamma_current * I_rV
	dD    <- delta * (I_wt + I_r + I_rV)
	dV    <- mu * RV + theta * S - beta * V + (I_r + I_rV) - w * V
	return(list(c(dS, dI_wt, dI_r, dI_rV, dR, dRV, dD, dV)))
}
