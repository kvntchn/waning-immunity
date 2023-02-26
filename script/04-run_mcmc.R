# mcmc.R
# February 21, 2023
# Fitting model via MCMC

library(R2jags)
san_francisco.dat <- fread("data/san_francisco.csv")
san_francisco.dat <- san_francisco.dat[date >= as.Date('2021-06-15')]
san_francisco.dat[,time := 1:.N]

source("script/mcmc.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN Metropolis-Hastings MCMC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(222)
mcmc_trace <- mh_mcmc(
	posterior = log_post_wrapper,
	init = parameters,
	constant_which = 1:6,
	num_iter = 6e4,
	progress = T,
	acceptance_progress = F)
# save(mcmc_trace, file = 'output/mcmc_trace.rdata')
plot(mcmc(mcmc_trace[,-(1:2)][-(1:2e4), -(1:6)]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ODE using parameter estimates from MCMC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load('output/mcmc_trace.rdata')
mcmc_trace <- as.data.table(mcmc_trace[,-(1:2)][-(1:3e4),])
mcmc_parameters <- colMeans(mcmc_trace[,-(1:2)])

# Initial parameters vs mcmc fit
data.frame(
	parameters = parameters,
	fitted = mcmc_parameters
)

times <- seq(from = 1, to = 603, by = parameters['dt'])
R0 <<- 0
trajectory.ode <- as.data.frame(ode(y = initial_state,
																		times = times,
																		parms = mcmc_parameters,
																		func = compartmental_model,
																		method = "lsode"))
rm(list = c('n_to_r', 'n_to_wt', 'R0'))

beta <- get_beta(trajectory.ode, parameters)
trajectory.ode$incidence <- beta * with(
	trajectory.ode, S * (I_wt + I_r + I_rV) + V * (I_r + I_rV))

# The first few entries of the trajectory matrix:
trajectory.ode %>%
	ggplot(aes(x = time)) +
	geom_point(aes(y = incidence), shape = 2, alpha = 0.6, col = 'salmon') +
	# geom_point(aes(y = I_wt + I_r + I_rV), shape = 2) +
	geom_point(data = san_francisco.dat[seq(1, .N, 3)],
						 aes(y = cases), size = 1/.pt) +
	coord_cartesian(ylim = c(0, 2e3)) +
	theme_bw()
