# ode.R
# February 9, 2023
# Simple ODE implementation

library(data.table)
library(deSolve)

san_francisco.dat <- fread("data/san_francisco.csv")


# Derivatives function for closed compartmental model:
compartmental_model <- function(time, state, parameters) {
	# Parameters:
	beta  <- parameters["R0"] / parameters["infectious_period"]
	gamma <- 1 / parameters["infectious_period"]
	delta <- parameters["death_rate"]
	theta0 <- parameters["vaccination_rate"]
	mu    <- 1 / parameters["recovery_period"]
	p     <- parameters['emergence_probability']
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
	# Time-varying parameter
	theta <- (1 - 0.25 / (S + R + I_wt + I_r)) * theta0 / (S + R + 0.05)
	# Derivatives:
	dS    <- mu * R - theta * S - beta / N * (I_wt + I_r + I_rV) * S
	dI_wt <- - (gamma + delta) * I_wt + beta / N * S * (I_r + I_rV)
	dI_r  <- - (gamma + delta) * I_r + beta / N * S * (I_r + I_rV)
	dI_rV <- - (gamma + delta) * I_rV + beta / N * V * (I_r + I_rV)
	dR    <- - mu * R - theta * R + gamma * (I_wt + I_r)
	dRV   <- - mu * RV + theta * R + gamma * I_rV
	dD    <- delta * (I_wt + I_r + I_rV)
	dV    <- mu * RV + theta * S - beta / N * V + (I_r + I_rV)
	return(list(c(dS, dI_wt, dI_r, dI_rV, dR, dRV, dD, dV)))
}

# Trajectory for model:
## Parameters
parameters <- c(
	R0 = 1.5,
	infectious_period = 14,
	death_rate = 1 / 1400,
	vaccination_rate = 1 / 100,
	recovery_period = 180)


# Initial conditions
N <- 8e5
initial_state <- c(
	S = N - 0.0001 * N,
	I_wt = 0.0001 * N,
	I_r = 0,
	R = 0,
	D = 0,
	V = 0,
	I_rV = 0,
	RV = 0)

# Times
times <- seq(from = 0, to = 930, by = 0.1)
# Produce ODE solution
trajectory <- ode(y = initial_state,
									times = times,
									parms = parameters,
									func = compartmental_model,
									method = "lsode")

# The first few entries of the trajectory matrix:
trajectory %>% as.data.frame() %>% ggplot(aes(x = time, y = I_wt + I_r + I_rV)) +
	geom_path() +
	theme_bw()
