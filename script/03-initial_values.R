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
	R0_low = 0.65,
	disease_duration = 14,
	recovery_period = 30 * 12,
	immune_period = 30 * 12,
	tighten_factor1 = 350,
	tighten_factor2 = 63,
	tighten_factor3 = 165,
	tighten_factor4 = 325,
	loosen_factor1 = 15,
	loosen_factor2 = 15,
	loosen_factor3 = 2.5,
	loosen_factor4 = 17
	# holiday_date =
	# 	# 160
	# 	as.integer(as.Date('2021-11-21') - as.Date('2021-06-14')),
	# reopening_date =
	# 	# 305
	# 	as.integer(as.Date('2022-04-15') - as.Date('2021-06-14')),
	# summer_date =
	# 	# 321
	# 	as.integer(as.Date('2022-05-01') - as.Date('2021-06-14'))
	)