# 05-simulations.R
# April 15, 2023
# Running model with heuristic choice of parameter values

library(lubridate)
library(foreach)
library(doParallel)

source('script/01-compartments.R')
san_francisco.dat <- fread("data/san_francisco.csv")
san_francisco.dat <- san_francisco.dat[date >= as.Date('2021-06-15')]
san_francisco.dat[,time := 1:.N]
source('script/03-initial_values.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation wrapper
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_sim <- function(
		parameters,
		times = 1:603,
		test_emergence_rate = NULL,
		test_immune_period = NULL) {
	current_frame <- sys.frame(sys.nframe())
	# cat('* '); print(current_frame)
	parameters['stochastic'] <- T
	parameters['emergence_date'] <- 0
	R0 <- parameters['R0_high1']
	resistant_strain_established <- F
	if (!is.null(test_emergence_rate)) {
		parameters['emergence_rate'] <- test_emergence_rate
	}
	if (!is.null(test_immune_period)) {
		parameters['immune_period'] <- test_immune_period
	}
	func <- function(time, state = initial_state, parameters) {
		compartmental_model(
			time, state, parameters, parent_frame = current_frame
		)}
	trajectory.ode <- as.data.frame(ode(
		y = initial_state,
		times = times,
		parms = parameters,
		func = func,
		method = "impAdams_d"))
	return(max(with(trajectory.ode, I_r + I_rV)))
}

# get_sim(parameters)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wrapper for simulation several times
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_sims <- function(
		B = 100, test_emergence_rate = NULL, test_immune_period = NULL,
		ncores = detectCores() - 2,
		param = parameters,
		one_sim = get_sim,
		quiet = T) {
	require(lubridate, quietly = T)
	require(foreach, quietly = T)
	require(doParallel, quietly = T)
	registerDoParallel(ncores)
	# start time
	start <- Sys.time()
	if (!quiet) {message("Run started at ", start, "\n")}
	sim <- foreach (b = 1:B) %dopar% {
		return(as.data.table(
			cbind(b = b,
						emergence_rate = test_emergence_rate,
						immune_period = test_immune_period,
						max_prev_resistant = one_sim(
							param,
							test_emergence_rate = test_emergence_rate,
							test_immune_period = test_immune_period))
		))
	}
	# end time
	time.unit <- "minutes"
	since_start <- lubridate::time_length(difftime(Sys.time(), start), time.unit)
	if (since_start > 90) {
		since_start <- since_start/60
		time.unit <- "hours"
	}
	if (!quiet) {
		cat(paste0("\n", round(since_start, 2), " ",
							 time.unit, " since get.bs() was called.\n"))
	}
	return(rbindlist(sim))
}

# get_sims()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run simulations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
emergence_rates <- seq(1e-4, 1e-2, length.out = 5)
immune_period <- c(1, 6, 12) * 30
param_grid <- expand.grid(
	emergence_rates = emergence_rates,
	immune_period = immune_period
)

# message(Sys.time())
# # 2023-04-21 14:38:57
# sim_output <- rbindlist(apply(param_grid, 1, function(param) {
# 	cat("Emergence rate:", param[[1]], "\n");
# 	cat("Immune period:", param[[2]], "\n");
# 	get_sims(1e3, test_emergence_rate = param[[1]], test_immune_period = param[[2]])
# }))
# message(Sys.time())
# # 2023-04-23 22:22:46
# # save(sim_output, file = 'output/sim_output.rdata')

# Establishment probabilities
sim_output[,.(establishment = mean(max_prev_resistant > 100) * 100), .(emergence_rate, immune_period)]
