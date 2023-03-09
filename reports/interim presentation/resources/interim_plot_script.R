# interim_plot_script
# March 6, 2023
# Running model with heuristic choice of parameter values; varying certain parameters and generating plots
# for interim presentation

library(tidyverse)
library(data.table)
library(deSolve)
library(gridExtra)

# Simple ODE implementation from
# source('~/Documents/capstone/capstone-epi/script/01-compartments.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Derivatives function for closed compartmental model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compartmental_model <- function(time, state = initial_state, parameters) {
  # Parameters:
  F_h     <- parameters['F_h']
  gamma   <- 1 / parameters["disease_duration"]
  delta   <- parameters["death_rate"]
  theta0  <- parameters["vaccination_rate"]
  w       <- 1 / parameters["immune_period"]
  mu      <- 1 / parameters["recovery_period"]
  p       <- 0
  R0_high <- parameters["R0_high"]
  R0_low  <- parameters["R0_low"]
  tighten_factor1   <- parameters["tighten_factor1"]
  tighten_factor2   <- parameters["tighten_factor2"]
  tighten_factor3   <- parameters["tighten_factor3"]
  tighten_factor4   <- parameters["tighten_factor4"]
  loosen_factor1   <- parameters["loosen_factor1"]
  loosen_factor2   <- parameters["loosen_factor2"]
  loosen_factor3   <- parameters["loosen_factor3"]
  loosen_factor4   <- parameters["loosen_factor4"]
  # as.integer(as.Date('2021-11-21') - as.Date('2021-06-14'))
  holiday_date     <- parameters["holiday_date"]
  # as.integer(as.Date('2022-04-15') - as.Date('2021-06-14'))
  reopening_date   <- parameters["reopening_date"]
  # as.integer(as.Date('2022-07-04') - as.Date('2021-06-14'))
  summer_date      <- parameters["summer_date"]
  emergence_date      <- parameters["emergence_date"]
  h       <- parameters["p_non_vax"]
  k       <- parameters["saturation"]
  dt      <- parameters["dt"]
  R0     <- get('R0', envir = .GlobalEnv)
  # States:
  S    <- state["S"]
  I_wt <- state["I_wt"]
  I_r  <- state['I_r']
  I_rV <- state["I_rV"]
  R    <- state["R"]
  RV   <- state["RV"]
  D    <- state["D"]
  V    <- state["V"]
  N <-  S + I_wt + I_r + I_rV + R + RV + D + V
  # Time-varying parameters
  ## Vaccination rate
  theta <- (1 - h / (S + R + I_wt + I_r)) * theta0 / (S + R + k)
  ## Time-varying force of  infection
  if (time > 0)   {F_h <- N / tighten_factor1; loosen_factor <- loosen_factor1}
  if (time >= holiday_date) {
    F_h <- N / tighten_factor2; loosen_factor <- loosen_factor2}
  if (time >= reopening_date) {
    F_h <- N / tighten_factor3; loosen_factor <- loosen_factor3}
  if (time >= summer_date) {
    F_h <- N / tighten_factor4; loosen_factor <- loosen_factor4}
  if ((I_wt + I_r + I_rV) > F_h & R0 > R0_low)                   {R0 <<- R0_low}
  if ((I_wt + I_r + I_rV) < (F_h / loosen_factor) & R0 < R0_high)  {R0 <<- R0_high}
  beta  <- get('R0', envir = .GlobalEnv) * gamma / N
  # Resistant strain?
  if ((I_r + I_rV) <= 1000 * N & time >= emergence_date) {
    p <- parameters['emergence_rate']
    if (parameters['stochastic']) {
      n_to_r  <<- rpois(1, dt * p * max(0, I_wt))
      n_to_wt <<- rpois(1, dt * p * max(0, I_r))
    } else {
      n_to_r  <<- dt * p * max(0, I_wt)
      n_to_wt <<- dt * p * max(0, I_r)
    }
    I_r  <- I_r  - n_to_wt + n_to_r
    I_wt <- I_wt + n_to_wt - n_to_r
  }
  # Derivatives:
  dS    <- mu * R + w * V - theta * S - beta * (I_wt + I_r + I_rV) * S
  dI_wt <- - (gamma + delta) * I_wt + beta * S * (I_wt)
  dI_r  <- (- (gamma + delta) * I_r  + beta * S * (I_r + I_rV))
  dI_rV <- - (gamma + delta) * I_rV + beta * V * (I_r + I_rV)
  dR    <- - mu * R - theta * R + gamma * (I_wt + I_r)
  dRV   <- - mu * RV + theta * R + gamma * I_rV
  dD    <- delta * (I_wt + I_r + I_rV)
  dV    <- mu * RV + theta * S - beta * V + (I_r + I_rV) - w * V
  return(list(c(dS, dI_wt, dI_r, dI_rV, dR, dRV, dD, dV)))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute what force of infection should be given state
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_beta <- function(traj, parameters) {
  beta <- rep(NA, nrow(traj))
  tighten_factor <- beta
  tighten_factor[traj[,'time'] > 0]    <- parameters['tighten_factor1']
  tighten_factor[traj[,'time'] >= 160] <- parameters['tighten_factor2']
  tighten_factor[traj[,'time'] >= 305] <- parameters['tighten_factor3']
  tighten_factor[traj[,'time'] >= 410] <- parameters['tighten_factor4']
  loosen_factor <- beta
  loosen_factor[traj[,'time'] > 0]    <- parameters['loosen_factor1']
  loosen_factor[traj[,'time'] >= 160] <- parameters['loosen_factor2']
  loosen_factor[traj[,'time'] >= 305] <- parameters['loosen_factor3']
  loosen_factor[traj[,'time'] >= 410] <- parameters['loosen_factor4']
  R0 <- c(beta, NA)
  R0[1] <- 0
  for (i in 1:length(beta)) {
    I <- (traj[i, 'I_wt'] + traj[i, 'I_r'] + traj[i, 'I_rV'])
    R0[i + 1] <- R0[i]
    if (I > N / tighten_factor[i] & R0[i] > parameters['R0_low']) {
      R0[i + 1] <- parameters['R0_low']}
    if (I < (N / tighten_factor[i] / loosen_factor[i]) & R0[i] < parameters['R0_high'])  {
      R0[i + 1] <- parameters['R0_high']}
  }
  return(R0[-1] / parameters['disease_duration'] / N)
}



# Starting values for ODE and parameters (chosen using epidemiologic knowledge) from
# source('~/Documents/capstone/capstone-epi/script/03-initial_values.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initial conditions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N <- 815201
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
  emergence_rate = 0,
  vaccination_rate = 250.6335,
  death_rate = 7.3e-4,
  dt = 1,
  stochastic = 0,
  holiday_date =  160.5, # as.integer(as.Date('2021-11-21') - as.Date('2021-06-14')),
  reopening_date =  305.5, # as.integer(as.Date('2022-04-15') - as.Date('2021-06-14')),
  summer_date =  385, # as.integer(as.Date('2022-07-04') - as.Date('2021-06-14'))
  emergence_date = 161.5, # as.integer(as.Date('2021-11-30') - as.Date('2021-06-14')),
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
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data Read-in
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

san_francisco.dat <- fread("~/Documents/capstone/data/san_francisco.csv")
san_francisco.dat <- san_francisco.dat[date >= as.Date('2021-06-15')]
san_francisco.dat[,time := 1:.N]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ODE solution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

times <- seq(from = 1, to = 603, by = parameters['dt'])
R0 <<- 0
trajectory.ode <- as.data.frame(ode(
  y = initial_state,
  times = times,
  parms = parameters,
  func = compartmental_model,
  method = "lsode"))
# rm(list = c('n_to_r', 'n_to_wt', 'R0'))

beta <- get_beta(trajectory.ode, parameters)
trajectory.ode$incidence <- beta * with(
  trajectory.ode, S * (I_wt + I_r + I_rV) + V * (I_r + I_rV))
trajectory.ode$incidence_wt <- beta * with(
  trajectory.ode, S * I_wt)
trajectory.ode$incidence_r <- beta * with(
  trajectory.ode, S * (I_r + I_rV) + V * (I_r + I_rV))

# The first few entries of the trajectory matrix:
#trajectory.ode %>%
#	ggplot(aes(x = time)) +
#	geom_path(aes(y = incidence, col = "Total"), shape = 2, alpha = 0.6) +
#	geom_path(aes(y = incidence_wt, col = "WT"), shape = 2, alpha = 0.6) +
#	geom_path(aes(y = incidence_r, col = "R"), shape = 2, alpha = 0.6) +
	# geom_path(aes(y = I_wt + I_r + I_rV), shape = 2) +
	# geom_point(data = san_francisco.dat,
	# 					 aes(y = cases), size = 1/.pt) +
#	geom_point(data = san_francisco.dat,
		#				 aes(y = deaths), size = 1/.pt) +
	# hospitalization rate?
#	coord_cartesian(ylim = c(0, 10)) +
#	theme_bw()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot from QMD file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The first few entries of the trajectory matrix:
#trajectory.ode %>%
#  ggplot(aes(x = time)) +
#  geom_path(aes(y = incidence), col = 'salmon') +
#  geom_point(
#    data = san_francisco.dat,
#    aes(y = cases), size = 1/.pt, alpha = 0.2) +
#  labs(subtitle = "R0_high = 2.2, R0_low = 0.65") +
#  theme_bw()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Increasing R_0 high and decreasing R_0 low
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Original Model: R0_high = 2.2, R0_low = 0.65

# function to change 'R0_high', 'R0_low'
pertb_R0 = function(R0_h, R0_l, init_params) {
  pertb_R0_params <- init_params
  pertb_R0_params[c('R0_high', 'R0_low')] <- c(R0_h, R0_l)
  R0 <<- 0
  set.seed(222)
  trajectory.odeR0 <- as.data.frame(ode(
    y = initial_state,
    times = times,
    parms = pertb_R0_params,
    func = compartmental_model,
    method = "lsode"))
  
  betaR0 <- get_beta(trajectory.odeR0, pertb_R0_params)
  trajectory.odeR0$incidence <- betaR0 * with(
    trajectory.odeR0, S * (I_wt + I_r + I_rV) + V * (I_r + I_rV))
  
  return(trajectory.odeR0)
}

# "in general, i would expect the rate of exponential growth to b larger for larger R0_high 
# and the rate of exponential decay to b faster for smaller R0_low"

# original covid strain R_0 = 2.79
# delta variant most dominant in CA june 2021, R_0 = 5.08; ranged from 3.2 to 8, with a mean of 5.08.
# BA. 4/5 strains have a R0 = 18.6
# The BQ.1.1 and BQ.1 variants continue to make up a majority of cases in Jan 2023, R_0 = 1.4 in Nov 2021


# increasing R0_h and fixing R0_l results in more rapid and frequent exponential growth; no prolonged waves; large, short-lived waves 
  # behavior is more predictable
plt11 <- ggplot(aes(x = time), data = trajectory.ode) +
  # initial model R0_h = 2.2, R0_l = 0.65
  geom_path(aes(y = incidence, col = 'Initial Model'), alpha = 0.65, lwd = 1.2, data = trajectory.ode) +
  # increasing R0_h and fixing R0_l 
  geom_path(aes(y = incidence, col = 'R_0 high = 3.0'), alpha = 0.5, lwd = 0.9,
            data = pertb_R0(R0_h = 3, R0_l = parameters['R0_low'], init_params = parameters)) + 
  # data
  geom_point(
    data = san_francisco.dat,
    aes(y = cases), size = 1/.pt, alpha = 0.2) +
  labs(title = "Increasing R_0 high, Fixing R_0 low") +
  scale_colour_discrete(name = "Model") + 
  theme_bw()

# fixing R0_h and decreasing R0_l results more frequent exponential growth during waves & major exponential growth in third wave
# R0_l = 0.3 (ranging from 0.15 to 0.53) : major exponential decay after exponential growth in third wave
# R0_l = 0.1 <' : prolonged third wave after initial exponential growth
plt12 <- ggplot(aes(x = time), data = trajectory.ode) +
  # initial model R0_h = 2.2, R0_l = 0.65
  geom_path(aes(y = incidence, col = 'Initial Model'), alpha = 0.65, lwd = 1, data = trajectory.ode) +
  # fix R0_h and decreasing R0_l to 0.3
  geom_path(aes(y = incidence, col = 'R_0 low = 0.3'), alpha = 0.4, lwd = 1.2,
            data = pertb_R0(R0_h = parameters['R0_high'], R0_l = 0.3, init_params = parameters)) + 
  # fix R0_h and decreasing R0_l to 0.1
  geom_path(aes(y = incidence, col = 'R_0 low = 0.1'), alpha = 0.65, lwd = 1.1,
            data = pertb_R0(R0_h = parameters['R0_high'], R0_l = 0.1, init_params = parameters)) + 
  # data
  geom_point(
    data = san_francisco.dat,
    aes(y = cases), size = 1/.pt, alpha = 0.2) +
  labs(title = "Fixing R_0 high, Decreasing R_0 low") +
  scale_colour_manual(name = "Model",
                      values = c('#F8766D', '#00BFC4', 'cyan4'),
                      breaks = c('Initial Model', 'R_0 low = 0.3', 'R_0 low = 0.1')) + 
  theme_bw()

grid.arrange(plt11, plt12, ncol=1)

ggsave(filename = 'plotR_0.png', 
       plot = grid.arrange(plt11, plt12, ncol=1),
       scale = 0.9,
       width = 7,
       height = 4,
       units = "in",
       device = "png",
       path = "/Users/laurendimaggio/Documents/capstone/plots")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Increasing Recovery Period/Duration of Immunity While Changing Tightening/Loosening Factors
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# duration of natural immunity = Recovery Period =  (days) (1/mu) = 1/(30*12)
# duration of vaccine protection = 'immune_period' (days) (1/omega) = 1/(30*12)

# function to change 'recovery_period', 'immune_period', 'tighten_factor', 'loosen_factor'
pertb_durat_tl_fact = function(pertb_params, init_params) {
  durat_tl_fact_params <- init_params
  param_selection <- tidyselect::contains(c('recovery_period', 'immune_period',
                                            'tighten_factor', 'loosen_factor'), 
                                          vars = names(durat_tl_fact_params))
  durat_tl_fact_params[param_selection] <- pertb_params
  R0 <<- 0
  set.seed(222)
  trajectory.ode <- as.data.frame(ode(
    y = initial_state,
    times = times,
    parms = durat_tl_fact_params,
    func = compartmental_model,
    method = "lsode"))
  
  beta <- get_beta(trajectory.ode, durat_tl_fact_params)
  trajectory.ode$incidence <- beta * with(
    trajectory.ode, S * (I_wt + I_r + I_rV) + V * (I_r + I_rV))
  
  return(trajectory.ode)
}

# increase immunity duration by 50% = immunity duration * 0.5
# decrease immunity duration by 50% = immunity duration * 1.5

# Increasing Recovery Period/Duration of Immunity + Increase Tightening/Loosening Factors
plt21 <- ggplot(aes(x = time), data = trajectory.ode) +
  # initial model 'recovery_period', 'immune_period', 'tighten_factor', 'loosen_factor'
  geom_path(aes(y = incidence, col = 'Initial Model'), alpha = 0.65, lwd = 1, data = trajectory.ode) +
  # Increasing Recovery Period/Duration of Immunity + Increase Tightening/Loosening Factors
    # increase recovery_period   immune_period to 540 days
  geom_path(aes(y = incidence, col = 'Increase mu, \nomega, h, l'), alpha = 0.55, lwd = 1.2,
            data = pertb_durat_tl_fact(pertb_params = c(parameters[c('recovery_period', 'immune_period')]*(12/18), 
                                                        375, 78, 190, 350, 25, 25, 4, 27), 
                                       init_params = parameters)) + 
  # data
  geom_point(
    data = san_francisco.dat,
    aes(y = cases), size = 1/.pt, alpha = 0.2) +
  labs(title = "Increase Tightening & Loosening Factors") + 
       #subtitle = "Increase Duration of Natural and Vaccinated Immunity by 2/3") +
  scale_colour_manual(name = "Model",
                      values = c('#F8766D', '#00BFC4'),
                      breaks = c('Initial Model', 'Increase mu, \nomega, h, l')
                      ) + 
  theme_bw()

# Increasing Recovery Period/Duration of Immunity + Decrease Tightening/Loosening Factors
plt22 <- ggplot(aes(x = time), data = trajectory.ode) +
  # initial model 'recovery_period', 'immune_period', 'tighten_factor', 'loosen_factor'
  geom_path(aes(y = incidence, col = 'Initial Model'), alpha = 0.7, lwd = 1, data = trajectory.ode) +
  # Increasing Recovery Period/Duration of Immunity + Decrease Tightening/Loosening Factors
    #  increase recovery_period   immune_period to 540 days
  geom_path(aes(y = incidence, col = 'Decrease h, l; \nincrease mu, \nomega'), alpha = 0.55, lwd = 1.2,
            data = pertb_durat_tl_fact(pertb_params = c(parameters[c('recovery_period', 'immune_period')]*(12/18),
                                                        325, 48, 140, 300, 9, 9, 1.5, 10), 
                                       init_params = parameters)) + 
  # data
  geom_point(
    data = san_francisco.dat,
    aes(y = cases), size = 1/.pt, alpha = 0.2) +
  labs(title = "Decrease Tightening & Loosening Factors") +
      # subtitle = "Increase Duration of Natural and Vaccinated Immunity by 2/3") +
  scale_colour_manual(name = "Model",
                      values = c('#F8766D', 'cyan3'),
                      breaks = c('Initial Model', 'Decrease h, l; \nincrease mu, \nomega')) + 
  theme_bw()

grid.arrange(plt21, plt22, ncol=1)

ggsave(filename = 'plot_TFL_F.png', 
       plot = grid.arrange(plt21, plt22, ncol=1),
       scale = 0.9,
       width = 7,
       height = 4,
       units = "in",
       device = "png",
       path = "/Users/laurendimaggio/Documents/capstone/plots")




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Varying Rates of Vaccination/Recovery/Immunity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# recovery_period = natural immunity (mu)
# immune_period = duration of vaccine protection (omega)
# vaccination_rate (theta)
# change 'vaccination_rate', 'recovery_period', 'immune_period'

# function to change 'vaccination_rate', 'recovery_period', 'immune_period'
pertb_vac_rec_immun = function(pertb_params, init_params) {
  vac_rec_immun_params <- init_params
  vac_rec_immun_params[c('vaccination_rate', 'recovery_period', 'immune_period')] <- pertb_params
  R0 <<- 0
  set.seed(222)
  trajectory.ode <- as.data.frame(ode(
    y = initial_state,
    times = times,
    parms = vac_rec_immun_params,
    func = compartmental_model,
    method = "lsode"))
  
  beta <- get_beta(trajectory.ode, vac_rec_immun_params)
  trajectory.ode$incidence <- beta * with(
    trajectory.ode, S * (I_wt + I_r + I_rV) + V * (I_r + I_rV))
  
  return(trajectory.ode)
}


# fixing 'recovery_period', 'immune_period', and increasing 'vaccination_rate'
plt3 <- ggplot(aes(x = time), data = trajectory.ode) +
  # initial model 'vaccination_rate', 'recovery_period', 'immune_period'
  geom_path(aes(y = incidence, col = 'Initial Model'), alpha = 0.8, lwd = 2, data = trajectory.ode) +
  # fixing 'recovery_period', 'immune_period', and increasing 'vaccination_rate' by 10% doesn't impact model much
  geom_path(aes(y = incidence, col = 'Increase theta by 10%; \nfix mu, omega'), alpha = 0.8, lwd = 0.7,
            data = pertb_vac_rec_immun(pertb_params = c(parameters['vaccination_rate']*1.1, 
                                                        parameters[c('recovery_period', 'immune_period')]), 
                                       init_params = parameters)) + 
  # fixing 'recovery_period', 'immune_period', and increasing 'vaccination_rate' by 50% impacts model during later wave
  geom_path(aes(y = incidence, col = 'Increase theta by 50%; \nfix mu, omega'), alpha = 0.7, lwd = 1,
            data = pertb_vac_rec_immun(pertb_params = c(parameters['vaccination_rate']*1.5, 
                                                        parameters[c('recovery_period', 'immune_period')]), 
                                       init_params = parameters)) + 
  # data
  geom_point(
    data = san_francisco.dat,
    aes(y = cases), size = 1/.pt, alpha = 0.2) +
  labs(title = "Fix Duration of Natural Immunity & Vaccine Protection; \nIncrease Vaccination Rate") +
  scale_colour_manual(name = "Model",
                      values = c('#F8766D', 'cyan4', 'cyan3'),
                      breaks = c('Initial Model', 'Increase theta by 10%; \nfix mu, omega', 
                                 'Increase theta by 50%; \nfix mu, omega')
                      ) + 
  theme_bw()

ggsave(filename = 'plot_incr_vax.png', 
       plot = plt3,
       scale = 0.9,
       width = 7,
       height = 4,
       units = "in",
       device = "png",
       path = "/Users/laurendimaggio/Documents/capstone/plots")

# (1) several distinct combinations of the parameters can lead to near identical fits
# (2) some choices lead to hugely different fits
# so on the one hand, we are under-identified for estimation, but on the other, we also have large sensitivity to changing some parameters
# so estimtion will be difficult 

##  due to reciprocals 
# increase immunity duration by 50% = immunity duration * 0.5
# decrease immunity duration by 50% = immunity duration * 1.5

# increasing theta by 50% is more powerful than halving immunity duration since vaccination was high to begin with

# vary 'recovery_period', 'immune_period',  increasing 'vaccination_rate'
plt31 <- ggplot(aes(x = time), data = trajectory.ode) +
  # initial model 'vaccination_rate', 'recovery_period', 'immune_period'
  geom_path(aes(y = incidence, col = 'Initial Model'), alpha = 0.8, lwd = 3, data = trajectory.ode) +
  # increase 'recovery_period', 'immune_period' by 50%, and increasing 'vaccination_rate' by 50% 
  geom_path(aes(y = incidence, col = 'Increase theta, mu, \nand omega by 50%'), alpha = 0.9, lwd = 2.2,
            data = pertb_vac_rec_immun(pertb_params = c(parameters['vaccination_rate']*1.5, 
                                                        parameters[c('recovery_period', 'immune_period')]*0.5), 
                                       init_params = parameters)) + 
  # fix 'recovery_period', 'immune_period'; increase 'vaccination_rate' by 50%
  geom_path(aes(y = incidence, col = 'Increase theta by 50%; \nfix mu, omega'), alpha = 0.8, lwd = 0.8,
            data = pertb_vac_rec_immun(pertb_params = c(parameters['vaccination_rate']*1.5, 
                                                        parameters[c('recovery_period', 'immune_period')]), 
                                       init_params = parameters)) + 
  # decrease 'recovery_period', 'immune_period' by 50%, and increasing 'vaccination_rate' by 50% 
  geom_path(aes(y = incidence, col = 'Increase theta by 50%; \ndecrease mu and \nomega by 50%'), alpha = 0.9, lwd = 0.6,
            data = pertb_vac_rec_immun(pertb_params = c(parameters['vaccination_rate']*1.5, 
                                                        parameters[c('recovery_period', 'immune_period')]*1.5), 
                                       init_params = parameters)) + 
  # decrease 'recovery_period', 'immune_period', 'vaccination_rate' by 50% 
  geom_path(aes(y = incidence, col = 'Decrease theta, mu, \nomega by 50%'), alpha = 0.9, lwd = 0.6,
            data = pertb_vac_rec_immun(pertb_params = c(parameters['vaccination_rate']*0.5, 
                                                        parameters[c('recovery_period', 'immune_period')]*1.5), 
                                       init_params = parameters)) + 
  # data
  geom_point(
    data = san_francisco.dat,
    aes(y = cases), size = 1/.pt, alpha = 0.2) +
  labs(title = "Vary Duration of Natural Immunity & Vaccine \nProtection; Increase Vaccination Rate by 50%") +
  scale_colour_manual(name = "Model",
                      #values = c('#F8766D', 'cyan3', 'cyan4'),
                      values = c('#F8766D', 'cyan2', '#E68613', 'cyan4', 'green'),
                      breaks = c('Initial Model', 'Increase theta, mu, \nand omega by 50%',
                                 'Increase theta by 50%; \nfix mu, omega',
                                 'Increase theta by 50%; \ndecrease mu and \nomega by 50%',
                                 'Decrease theta, mu, \nomega by 50%')
                      ) + 
  theme_bw()

# varying 'vaccination_rate', and let 'recovery_period' = 3 mos, 'immune_period' = 5 mos 
plt32 <- ggplot(aes(x = time), data = trajectory.ode) +
  # initial model 'vaccination_rate', 'recovery_period', 'immune_period'
  geom_path(aes(y = incidence, col = 'Initial Model'), alpha = 0.8, lwd = 2, data = trajectory.ode) +
  # fixing 'vaccination_rate', and let 'recovery_period' = 3 mos, 'immune_period' = 5 mos doesn't impact model much
  geom_path(aes(y = incidence, col = 'Fix theta'), alpha = 0.7, lwd = 1.3,
            data = pertb_vac_rec_immun(pertb_params = c(parameters['vaccination_rate'], 30*3, 30*50), 
                                       init_params = parameters)) + 
  # decrease 'vaccination_rate' by 50%, let 'recovery_period' = 3 mos, 'immune_period' = 5 mos impacts model some during later waves
  geom_path(aes(y = incidence, col = 'Decrease theta \nby 50%'), alpha = 0.9, lwd = 0.8,
            data = pertb_vac_rec_immun(pertb_params = c(parameters['vaccination_rate']*0.5, 30*3, 30*50), 
                                       init_params = parameters)) + 
  # increase 'vaccination_rate' by 50% (or +70%), let 'recovery_period' = 3 mos, 'immune_period' = 5 mos dampens third wave
  geom_path(aes(y = incidence, col = 'Increase theta \nby 50%'), alpha = 0.9, lwd = 0.7,
            data = pertb_vac_rec_immun(pertb_params = c(parameters['vaccination_rate']*1.5, 30*3, 30*50), 
                                       init_params = parameters)) + 
  # data
  geom_point(
    data = san_francisco.dat,
    aes(y = cases), size = 1/.pt, alpha = 0.2) +
  labs(title = "Vary Vaccination Rate; Fix Duration of Natural \nImmunity and Vaccine Protection at 3 Months \nand 5 Months") +
  scale_colour_manual(name = "Model",
                      values = c('#F8766D', '#E68613', 'cyan2', 'cyan4'),
                      breaks = c('Initial Model', 'Fix theta', 'Decrease theta \nby 50%', 'Increase theta \nby 50%')) + 
  theme_bw()


grid.arrange(plt31, plt32, ncol=2)

ggsave(filename = 'plot_vac_immun.png', 
       plot = grid.arrange(plt31, plt32, ncol=2),
       scale = 7/6,
       width = 8,
       height = 4,
       units = "in",
       device = "png",
       path = "/Users/laurendimaggio/Documents/capstone/plots")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Varying Rates of Emergence
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###  emergence rate p: 1 in a million should b small, but i think 1 in a thousand would b large

pertb_p = function(pertb_param, init_params) {
  pert_p_param <- init_params
  pert_p_param['emergence_rate'] <- pertb_param
  R0 <<- 0
  set.seed(222)
  trajectory.ode <- as.data.frame(ode(
    y = initial_state,
    times = times,
    parms = pert_p_param,
    func = compartmental_model,
    method = "lsode"))
  
  beta <- get_beta(trajectory.ode, pert_p_param)
  trajectory.ode$incidence <- beta * with(
    trajectory.ode, S * (I_wt + I_r + I_rV) + V * (I_r + I_rV))
  
  return(trajectory.ode)
}

# small p = 1/million doesn't seem to impact model 
plt41 <- ggplot(aes(x = time), data = trajectory.ode) +
  # initial model p = 0
  geom_path(aes(y = incidence, col = 'Initial Model'), alpha = 0.9, lwd = 2.5, data = trajectory.ode) +
  # small p = 1/million
  geom_path(aes(y = incidence, col = 'p = 1/1e6'), alpha = 0.8, lwd = 0.6,
            data = pertb_p(pertb_param = 1/1e6, init_params = parameters)) + 
  # data
  geom_point(data = san_francisco.dat, aes(y = cases), size = 1/.pt, alpha = 0.2) +
  labs(title = "Small p") +
  scale_colour_manual(name = "Model",
                      values = c('#F8766D', 'cyan2')) +
  theme_bw()

# large p = 1/thousand doesn't barely impacts model 
plt42 <- ggplot(aes(x = time), data = trajectory.ode) +
  # initial model p = 0
  geom_path(aes(y = incidence, col = 'Initial Model'), alpha = 0.9, lwd = 2.5, data = trajectory.ode) +
  # large p = 1/thousand
  geom_path(aes(y = incidence, col = 'p = 1/1000'), alpha = 0.8, lwd = 0.6,
            data = pertb_p(pertb_param = 1/1000, init_params = parameters)) + 
  # data
  geom_point(data = san_francisco.dat, aes(y = cases), size = 1/.pt, alpha = 0.2) +
  labs(title = "Large p") +
  scale_colour_manual(name = "Model",
                      values = c('#F8766D', 'cyan2')) +
  theme_bw()

grid.arrange(plt41, plt42, ncol=2)

ggsave(filename = 'plot_p.png', 
       plot = grid.arrange(plt41, plt42, ncol=2),
       scale = 5/4,
       width = 8,
       height = 4,
       units = "in",
       device = "png",
       path = "/Users/laurendimaggio/Documents/capstone/plots")


