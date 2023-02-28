# mcmc.R
# February 21, 2023
# Fitting model via MCMC

library(data.table)
library(R2jags)
san_francisco.dat <- fread("data/san_francisco.csv")
san_francisco.dat <- san_francisco.dat[date >= as.Date('2021-06-15')]
san_francisco.dat[,time := 1:.N]

source("script/02-mcmc.R")
source('script/03-initial_values.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN Metropolis-Hastings MCMC (No emergence)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prior_mean_no_emergence <- parameters[-(1:11)]
prior_sd_no_emergence   <- c(
	0.2, 0.1,
	14,
	30 * 6, 30 * 6,
	500, 100, 500, 500,
	5, 5, 2, 5)
prior_hyperparameters <- data.frame(
	param = names(prior_mean_no_emergence),
	alpha = prior_mean_no_emergence^2 / prior_sd_no_emergence^2,
	beta  = prior_mean_no_emergence / prior_sd_no_emergence^2
)

set.seed(222)
mcmc_trace_no_emergence <- mh_mcmc(
	posterior = log_post_wrapper,
	init = parameters,
	constant_which = 1:11,
	num_iter = 6e4,
	progress = T,
	acceptance_progress = F)
# save(mcmc_trace_no_emergence, file = 'output/mcmc_trace_no_emergence.rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ODE using parameter estimates from MCMC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load('output/mcmc_trace_no_emergence.rdata')
# plot(mcmc(mcmc_trace_no_emergence[,-(1:2)][, -(1:11)]))
# plot(mcmc(mcmc_trace_no_emergence[,-(1:2)][-(1:2e4), -(1:11)]))
mcmc_trace <- as.data.table(mcmc_trace_no_emergence[,-(1:2)][-(1:3e4),])
mcmc_parameters <- colMeans(mcmc_trace)

# Initial parameters vs mcmc fit
data.frame(
	prior_mean = parameters,
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
