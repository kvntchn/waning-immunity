# priors.R
# February 21, 2023
# Priors for state model parameters

source('script/compartments.R')
san_francisco.dat <- fread("data/san_francisco.csv")
san_francisco.dat <- san_francisco.dat[-(1:75),]
san_francisco.dat[,time := 1:.N]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initial conditions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parameters <- c(
	F_h = N / 2000,
	R0_high = 1.5,
	R0_low = 0.75,
	R0_omicron = 3,
	R0_previous = -Inf,
	immune_period = Inf,
	disease_duration = 17,
	disease_duration_omicron = 15,
	death_rate = 1 / 1000,
	vaccination_rate = 1 / 1000,
	recovery_period = 30 * 12,
	recovery_period_omicron = 30 * 10,
	p_non_vax = 0.2,
	saturation = 0.01,
	emergence_probability = 5e-4,
	omicron_date = 390,
	dt = 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ODE solution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
times <- seq(from = 0, to = 930, by = parameters['dt'])
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

