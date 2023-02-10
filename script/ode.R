# ode.R
# February 9, 2023
# Simple ODE implementation

library(tidyverse)
library(data.table)
library(deSolve)

san_francisco.dat <- fread("data/san_francisco.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Derivatives function for closed compartmental model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compartmental_model <- function(time, state, parameters) {
	# Parameters:
	F_h     <- parameters['SIP']
	gamma   <- 1 / parameters["disease_duration"]
	delta   <- parameters["death_rate"]
	theta0  <- parameters["vaccination_rate"]
	mu      <- 1 / parameters["recovery_period"]
	p       <- parameters['emergence_probability']
	R0_high <- parameters["R0_high"]
	R0_low  <- parameters["R0_low"]
	if (!exists('R0_prev')) {R0_prev <- parameters["R0_previous"]}
	h       <- parameters["p_non_vax"]
	k       <- parameters["saturation"]
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
	## Vaccination rate
	theta <- (1 - h / (S + R + I_wt + I_r)) * theta0 / (S + R + k)
	## Time-varying force of  infection
	if (I > F_h   & R0_prev > R0_low)  {R0 <<- R0_low}
	if (I < F_h/8 & R0_prev < R0_high) {R0 <<- R0_high}
	beta  <- R0 * gamma / N
	R0_prev <<- R0
	# Derivatives:
	dS    <- mu * R - theta * S - beta * (I_wt + I_r + I_rV) * S
	dI_wt <- - (gamma + delta) * I_wt + beta * S * (I_wt)
	dI_r  <- - (gamma + delta) * I_r  + beta * S * (I_r + I_rV)
	dI_rV <- - (gamma + delta) * I_rV + beta * V * (I_r + I_rV)
	dR    <- - mu * R - theta * R + gamma * (I_wt + I_r)
	dRV   <- - mu * RV + theta * R + gamma * I_rV
	dD    <- delta * (I_wt + I_r + I_rV)
	dV    <- mu * RV + theta * S - beta * V + (I_r + I_rV)
	return(list(c(dS, dI_wt, dI_r, dI_rV, dR, dRV, dD, dV)))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Trajectory for model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initial conditions
N <- 8e5
initial_state <- c(
	S = N - 100,
	I_wt = 100,
	I_r = 0,
	I_rV = 0,
	R = 0,
	RV = 0,
	D = 0,
	V = 0)

# Parameters
parameters <- c(
	SIP = N / 500,
	R0_high = 2.5,
	R0_low = 0.77,
	R0_previous = -Inf,
	disease_duration = 14/0.99,
	death_rate = 1 / 1000,
	vaccination_rate = 1 / 1000,
	recovery_period = 180,
	p_non_vax = 0.01,
	saturation = 0.01)


# Times
times <- seq(from = 0, to = 930, by = 0.1)
# Produce ODE solution
R0_log <- c()
trajectory.ode <- ode(y = initial_state,
											times = times,
											parms = parameters,
											func = compartmental_model,
											method = "lsode")

# The first few entries of the trajectory matrix:
trajectory.ode %>% as.data.frame() %>%
	ggplot(aes(x = time)) +
	geom_path(aes(y = I_wt), col = 'forestgreen') +
	geom_path(aes(y = I_r + I_rV), col = 'salmon') +
	theme_bw()

