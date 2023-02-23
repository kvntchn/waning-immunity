# eyeball-fit.R
# February 21, 2023
# Running model with heuristic choice of parameter values

source('script/compartments.R')
source('script/initial_values.R')
san_francisco.dat <- fread("data/san_francisco.csv")
san_francisco.dat <- san_francisco.dat[date >= as.Date('2021-06-15')]
san_francisco.dat[,time := 1:.N]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ODE solution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
times <- seq(from = 1, to = 603, by = parameters['dt'])
R0 <<- 0
set.seed(222)
trajectory.ode <- as.data.frame(ode(y = initial_state,
											times = times,
											parms = parameters,
											func = compartmental_model,
											method = "lsode"))
rm(list = c('n_to_r', 'n_to_wt', 'R0'))

beta <- get_beta(trajectory.ode, parameters)
trajectory.ode$incidence <- beta * with(
	trajectory.ode, S * (I_wt + I_r + I_rV) + V * (I_r + I_rV))

# The first few entries of the trajectory matrix:
trajectory.ode %>%
	ggplot(aes(x = time)) +
	geom_point(aes(y = incidence), shape = 2) +
	# geom_point(aes(y = I_wt + I_r + I_rV), shape = 2) +
	geom_point(data = san_francisco.dat,
						aes(y = cases), size = 1/.pt) +
	coord_cartesian(ylim = c(0, 2e3)) +
	theme_bw()