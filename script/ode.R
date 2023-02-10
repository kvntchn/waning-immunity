# ode.R
# February 9, 2023
# Simple ODE implementation

library(tidyverse)
library(data.table)
library(deSolve)

san_francisco.dat <- fread("data/san_francisco.csv")
san_francisco.dat <- san_francisco.dat[-(1:75),]
san_francisco.dat[,time := 1:.N]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Derivatives function for closed compartmental model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compartmental_model <- function(time, state, parameters) {
	# Parameters:
	F_h     <- parameters['F_h']
	gamma   <- 1 / parameters["disease_duration"]
	delta   <- parameters["death_rate"]
	theta0  <- parameters["vaccination_rate"]
	mu      <- 1 / parameters["recovery_period"]
	p       <- parameters['emergence_probability']
	R0_high <- parameters["R0_high"]
	R0_low  <- parameters["R0_low"]
	R0_holiday  <- parameters["R0_holiday"]
	if (!exists('R0_prev')) {R0_prev <- parameters["R0_previous"]}
	h       <- parameters["p_non_vax"]
	k       <- parameters["saturation"]
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
	## Vaccination rate
	theta <- (1 - h / (S + R + I_wt + I_r)) * theta0 / (S + R + k)
	## Time-varying force of  infection
	if (I > F_h   & R0_prev > R0_low)  {R0 <<- R0_low}
	if (I < F_h/8 & R0_prev < R0_high) {R0 <<- R0_high}
	beta  <- R0 * gamma / N
	R0_prev <<- R0
	# Resistant strain?
	if ((I_r + I_rV) <= 0.0001 * N) {
		n_to_r  <<- rpois(1, dt * p * max(0, I_wt))
		n_to_wt <<- rpois(1, dt * p * max(0, I_r))
	}
	I_r  <- I_r  - n_to_wt + n_to_r
	I_wt <- I_wt + n_to_wt - n_to_r
	# California goes wild
	if (time >= 409 & time <= 444) {R0 <<- R0_holiday}
	if (time >= 555 & time <= 572) {R0 <<- R0_holiday}
	# Derivatives:
	dS    <- mu * R - theta * S - beta * (I_wt + I_r + I_rV) * S
	dI_wt <- - (gamma + delta) * I_wt + beta * S * (I_wt)
	dI_r  <- - (gamma + delta) * I_r  + beta * S * (I_r + I_rV) + p * I_wt
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
	S = N - 30,
	I_wt = 30,
	I_r = 0,
	I_rV = 0,
	R = 0,
	RV = 0,
	D = 0,
	V = 0)

# Parameters
parameters <- c(
	F_h = N / 2000,
	R0_high = 1.5,
	R0_low = 0.75,
	R0_holiday = 3,
	R0_previous = -Inf,
	disease_duration = 15,
	death_rate = 1 / 1000,
	vaccination_rate = 1 / 1000,
	recovery_period = 30 * 12,
	p_non_vax = 0.2,
	saturation = 0.01,
	emergence_probability = 2e-4,
	dt = 1)


# Times
times <- seq(from = 0, to = 930, by = parameters['dt'])
# Produce ODE solution
R0_log <- c()
set.seed(222)
trajectory.ode <- ode(y = initial_state,
											times = times,
											parms = parameters,
											func = compartmental_model,
											method = "lsode")

# The first few entries of the trajectory matrix:
trajectory.ode %>% as.data.frame() %>%
	ggplot(aes(x = time)) +
	geom_path(aes(y = I_wt + I_r + I_rV), alpha = 0.3) +
	geom_path(aes(y = I_wt), col = 'forestgreen') +
	geom_path(aes(y = I_r + I_rV), col = 'salmon') +
	geom_point(data = san_francisco.dat,
						aes(y = cases), size = 1/.pt) +
	theme_bw()

