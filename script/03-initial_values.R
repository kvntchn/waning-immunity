# initial_values.R
# February 22, 2023
# Starting values for ODE and parameters (chosen using epidemiologic knowledge)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initial conditions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N <- 815201
initial_state <- c(
	S = 815201 - 100,
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

# Time-varying
vaccination_rate.loess <- loess(fully_vaccinated ~ time, san_francisco.dat, span = 0.6)
# plot(san_francisco.dat$time, predict(vaccination_rate.loess))
# points(san_francisco.dat$time, san_francisco.dat$fully_vaccinated)
vaccination_rate <- data.frame(
	times = san_francisco.dat$time,
	theta = shift(predict(vaccination_rate.loess), 14, fill = 0) / N
)
vaccination_rate_input <- approxfun(vaccination_rate, rule = 2)

death_rate.loess <- loess(deaths ~ time, san_francisco.dat, span = 0.6)
# plot(san_francisco.dat$time, predict(death_rate.loess))
death_rate <- data.frame(
	times = san_francisco.dat$time,
	delta = predict(death_rate.loess) / N
)
death_rate_input <- approxfun(death_rate, rule = 2)

parameters <- c(
	p_non_vax = 0.1853923,
	saturation = 0.01,
	emergence_rate = 0,
	dt = 1,
	stochastic = 0,
	emergence_date = 161.5, # as.integer(as.Date('2021-11-30') - as.Date('2021-06-14')),
	# Random parameters
	holiday_date =  176, # as.integer(as.Date('2021-11-23') - as.Date('2021-06-14')),
	reopening_date =  298, # as.integer(as.Date('2022-04-15') - as.Date('2021-06-14')),
	summer_date =  399, # as.integer(as.Date('2022-07-04') - as.Date('2021-06-14'))
	R0_high1 = 2.1,
	R0_high2 = 2.1,
	R0_high3 = 1.5,
	R0_high4 = 1.7,
	R0_low1 = 0.75,
	R0_low2 = 0.5,
	R0_low3 = 0.75,
	R0_low4 = 0.75,
	disease_duration = 14,
	recovery_period = 30 * 12,
	immune_period = 30 * 12,
	tighten_factor1 = 250,
	tighten_factor2 = 40,
	tighten_factor3 = 150,
	tighten_factor4 = 350,
	loosen_factor1 = 12,
	loosen_factor2 = 15,
	loosen_factor3 = 1.2,
	loosen_factor4 = 5,
	NULL
)

# parameters[-constant_which] <-
# 	c(
# 		173.5918,
# 		304.3872,
# 		401.7339,
# 		2.086288,
# 		2.059187,
# 		1.420506,
# 		1.669773,
# 		0.7466942,
# 		0.5209631,
# 		0.7483176,
# 		0.7419277,
# 		13.99949,
# 		309.8113,
# 		33.24851,
# 		155.972,
# 		35.58448,
# 		99.05499,
# 		295.4496,
# 		13.27637,
# 		15.75084,
# 		1.412681,
# 		4.323944
# 	)
