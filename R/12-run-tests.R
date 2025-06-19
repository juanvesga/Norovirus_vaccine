# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root           <- here::here()
infile0        <- file.path(root,"output", "parameters.qs2")
infile_dataS   <- file.path(root,"output", "data_short.qs2")
infile_dataL   <- file.path(root,"output", "data_long.qs2")
infile_input   <- file.path(root,"output", "params_list.qs2")
infile_prior   <- file.path(root,"output", "priors.qs2")
infile_model2  <- file.path(root,"models", "model_2pars.R")
infile_chains  <- file.path(root, "output", "chains.qs")
infile_index   <- file.path(root, "output", "index.qs")
infile_odin1chain   <- file.path(root, "output", "chains_imm2par.RData")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))

# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))

library(odin2)
library(dust2)
library(monty)
library(posterior)
library(bayesplot)

#
#
#

update_parameters<-function(ins,params){
  
  params$beta_1     = ins["beta_1"]
  params$beta_2     = ins["beta_2"]
  params$beta_3     = ins["beta_3"]
  params$beta_4     = ins["beta_4"]
  params$aduRR      = ins["aduRR"]
  params$maternalAB = ins["maternalAB"]
  params$imm_yr     = ins["imm_yr"]
  params$imm_fac    = ins["imm_fac"]
  params$repfac_0   = ins["repfac_0"]
  params$repfac_5   = ins["repfac_5"]
  params$repfac_15  = ins["repfac_15"]
  params$repfac_65p = ins["repfac_65p"]
  params$crossp_GI  = ins["crossp_GI"]
  params$crossp_GII = ins["crossp_GII"]

return(params)
}

# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------


# Source Odin model
source(infile_model2)
parameters <- qs_read(infile0)
data       <- qs_read(infile_dataS)
pars_list  <- qs_read(infile_input)
prior      <- qs_read(infile_prior)
#chains     <- qs_read(infile_chains) 
chains_odin1<-load(infile_odin1chain)
calibrated_pars<-exp(processed_chains$pars)
init_states    <-processed_chains$state
init_states    <-init_states[-c(1),]  ## Remove time vector from initial state

rm(processed_chains)


# Set staring parameters
pars<-update_parameters(calibrated_pars[1,],pars_list$pars_list)  

#dust object
sys  <- dust_system_create(noro_model(), pars, deterministic = TRUE, n_particles = 1)

# Get state index
index<- dust_unpack_index(sys)

#Get default state initial
state<-dust_system_state(sys)

named_state<-dust_unpack_state(sys,state)

# Set inital model state based on calibrated model

# Set default initial state of the model 
dust_system_set_state_initial(sys)

#init0<-init_states[,1]

#names(init0)<-paste0(names(named_state))
  
#dust_system_set_state(sys,init0)

# Set parameters
dust_system_update_pars(sys, pars=pars)


# run simulation
tt<-seq(0,24109+365*10)
y0<- dust_system_simulate(sys, tt)

# unpack results
y0 <- dust_unpack_state(sys, y0)

plot(tt,(y0$new_cases),
        type="l",
        xlim = c(24109,24109+365*10),
        ylim = c(0, 1e5))


plot(tt[tt %% 7 == 0], y0$new_cases_week[tt %% 7 == 0], type = "l",
     ylab = "Weekly incidence", 
     xlab = "Time",
     xlim = c(24109,24109+365*10),
     ylim=c(0,5e5))







# Generate Initial conditions vector
seed<-1# Number infected per strain
#init<-generate_initial_state(index,state,parameters,seed,2)


# check initial conditions in current object
dust_system_state(sys)

# Pass initial conditions to model
dust_system_set_state_initial(sys)

# Check again initial conditions in object
dust_system_state(sys)

#Simulate system up to a point 
t <- seq(0, 20000)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)


matplot(t(y$pop_all),type='l')

plot(y$pop_all[8,])

plot(t, y$new_cases_week, type = "l", xlab = "Time", ylab = "Infected population",
     xlim = c(1e4,2e4),ylim=c(0,5e4))
# 
plot(t[t %% 7 == 0], y$new_cases_week[t %% 7 == 0], type = "o", pch = 19,
     ylab = "Infection incidence", xlab = "Time",xlim = c(1e4,10600),ylim=c(0,15e4))
lines(t, y$new_cases_week, col = "red")