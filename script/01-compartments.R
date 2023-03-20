# state_model.R
# February 9, 2023
# Simple ODE implementation

library(tidyverse)
library(data.table)
library(deSolve)
library(rlang)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Derivatives function for closed compartmental model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compartmental_model <- function(time, state = initial_state, parameters) {
	# Parameters:
	with(as.list(c(parameters, state)), {
		delta   <- death_rate_input(time)
		theta0  <- vaccination_rate_input(time)
		N <-  S + I_wt + I_r + I_rV + R + RV + D + V
		## Vaccination rate
		theta <- (1 - p_non_vax / (S + R + I_wt + I_r)) * theta0 / (S + R + saturation)

		## Time-varying force of  infection
		time_varying_parameters <- data.frame(
			dates = c(0, holiday_date, reopening_date, summer_date),
			tighten_factor = c(tighten_factor1, tighten_factor2,
												 tighten_factor3, tighten_factor4),
			loosen_factor = c(loosen_factor1, loosen_factor2,
												loosen_factor3, loosen_factor4),
			R0_high = c(R0_high1, R0_high2, R0_high3, R0_high4),
			R0_low = c(R0_low1, R0_low2, R0_low3, R0_low4)
		)
		tighten_factor <- approxfun(
			time_varying_parameters[,c(1, 2)], method = 'constant', rule = 2)
		loosen_factor <- approxfun(
			time_varying_parameters[,c(1, 3)], method = 'constant', rule = 2)
		R0_high <- approxfun(
			time_varying_parameters[,c(1, 4)], method = 'linear', rule = 2)
		R0_low <- approxfun(
			time_varying_parameters[,c(1, 5)], method = 'linear', rule = 2)

		env <- current_env()
		# if (all(sapply(c("I_wt", "I_r", "I_rV", "N", "tighten_factor", "loosen_factor"),
		# 							 function(x) {is.finite(get(x, envir = env))}))) {
			if ((I_wt + I_r + I_rV) > (N / tighten_factor(time)) &
					get('R0', envir = .GlobalEnv) > 1) {
				assign('R0', R0_low(time), envir = .GlobalEnv)
			} else if ((I_wt + I_r + I_rV) < (N / (tighten_factor(time) * loosen_factor(time))) &
								 get('R0', envir = .GlobalEnv) < 1) {
				assign('R0', R0_high(time), envir = .GlobalEnv)}

			# Resistant strain?
			if ((I_r + I_rV) >= N / 1000 | get("resistant_strain_established", .GlobalEnv)) {
				resistant_strain_established <<- T
			} else if (time >= emergence_date) {
				if (stochastic) {
					n_to_r  <<- rpois(1, dt * emergence_rate * max(0, I_wt))
					n_to_wt <<- rpois(1, dt * emergence_rate * max(0, I_r))
				} else {
					n_to_r  <<- dt * emergence_rate * max(0, I_wt)
					n_to_wt <<- dt * emergence_rate * max(0, I_r)
				}
				I_r  <- I_r  - n_to_wt + n_to_r
				I_wt <- I_wt + n_to_wt - n_to_r
			}
		# }

		beta  <- get('R0', envir = .GlobalEnv) / disease_duration / N

		# Derivatives:
		dS    <- (1 / recovery_period) * R +
			(1 / immune_period) * V - theta * S - beta * (I_wt + I_r + I_rV) * S
		dI_wt <- - ((1 / disease_duration) + delta) * I_wt  + beta * S * (I_wt)
		dI_r  <- (- ((1 / disease_duration) + delta) * I_r  + beta * S * (I_r + I_rV))
		dI_rV <- - ((1 / disease_duration) + delta) * I_rV + beta * V * (I_r + I_rV)
		dR    <- - (1 / recovery_period) * R - theta * R + (1 / disease_duration) * (I_wt + I_r)
		dRV   <- - (1 / recovery_period) * RV + theta * R + (1 / disease_duration) * I_rV
		dD    <- delta * (I_wt + I_r + I_rV)
		dV    <- (1 / recovery_period) * RV + theta * S - beta * V + (I_r + I_rV) -
			(1 / immune_period) * V
		return(list(
			c(dS, dI_wt, dI_r, dI_rV, dR, dRV, dD, dV),
			death_rate = delta,
			vaccination_rate = theta0,
			R0 = R0,
			beta = unname(beta)
			# tighten_factor = tighten_factor,
			# loosen_factor  = loosen_factor
		))
	})
}