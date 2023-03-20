# eyeball-fit.R
# February 21, 2023
# Running model with heuristic choice of parameter values

source('script/01-compartments.R')
san_francisco.dat <- fread("data/san_francisco.csv")
san_francisco.dat <- san_francisco.dat[date >= as.Date('2021-06-15')]
san_francisco.dat[,time := 1:.N]
source('script/03-initial_values.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ODE solution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
times <- seq(from = 1, to = nrow(san_francisco.dat), by = parameters['dt'])
R0 <- parameters['R0_high1']
resistant_strain_established <- F
trajectory.ode <- as.data.frame(ode(
	y = initial_state,
	times = times,
	parms = parameters,
	func = compartmental_model,
	method = "impAdams_d"))
rm(list = c('n_to_r', 'n_to_wt', 'R0'))

trajectory.ode$incidence_wt <- with(
	trajectory.ode, beta * S * I_wt)
trajectory.ode$incidence_r <- with(
	trajectory.ode, beta * (S * (I_r + I_rV) + V * (I_r + I_rV)))
trajectory.ode$incidence <- with(
	trajectory.ode, incidence_wt + incidence_r)

# The first few entries of the trajectory matrix:
trajectory.ode %>%
	ggplot(aes(x = time)) +
	geom_path(aes(y = incidence, col = "Total"), alpha = 0.6) +
	# geom_path(aes(y = incidence_wt, col = "WT"), alpha = 0.6) +
	# geom_path(aes(y = incidence_r, col = "R"), alpha = 0.6) +
	# geom_path(aes(y = I_wt + I_r + I_rV), alpha = 0.6) +
	geom_point(data = san_francisco.dat,
						 aes(y = cases), size = 1/.pt) +
	# geom_point(data = san_francisco.dat,
	# 					 aes(y = deaths), size = 1/.pt) +
	# hospitalization rate?
	coord_cartesian(ylim = c(0, max(san_francisco.dat$cases))) +
	theme_bw()
