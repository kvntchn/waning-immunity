# initial_values.R
# February 22, 2023
# Starting values for ODE and parameters (chosen using epidemiologic knowledge)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initial conditions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parameters <- c(
	p_non_vax = 0.2,
	saturation = 0.01,
	emergence_probability = 0,
	vaccination_rate = 250.6335,
	death_rate = 7.3e-4,
	dt = 1,
	# Random parameters
	R0_high = 2.2,
	R0_low = 0.7,
	disease_duration = 14,
	recovery_period = 30 * 12,
	immune_period = 30 * 12,
	tighten_factor1 = 300,
	tighten_factor2 = 65,
	tighten_factor3 = 120,
	tighten_factor4 = 220,
	loosen_factor1 = 15,
	loosen_factor2 = 12,
	loosen_factor3 = 2,
	loosen_factor4 = 12
	)