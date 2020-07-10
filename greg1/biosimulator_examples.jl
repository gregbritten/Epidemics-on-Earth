
using BioSimulator
using Plots


############################################################
### BIRTH DEATH ############################################
############################################################
network = Network("BD") #initialize network

##--DEFINE STATE VARIABLES AS 'SPECIES'--#####
network <= Species("X", 5)

##--DEFINE REACTIONS--###########
network <= Reaction("birth", 2.0, "X --> X + X")
network <= Reaction("death", 1.0, "X --> 0")

##--RUN SIMULATION--############
result = simulate(model, Direct(), tfinal = 5.0, rates_cache = HasRates, save_points = nothing, ntrials=5)

##--PLOT--#############
p = plot(result[1])
for i=2:length(result)
    p = plot!(result[i], summary = :trajectory)
end
display(p)


########################################################
### SIR ################################################
########################################################
n_people = 1000
β₀       = 0.34
β        = β₀/n_people
γ        = 0.25
days     = 100

SIR = Network("SIR")
SIR <= Species("S",n_people - 1)
SIR <= Species("I",1)
SIR <= Species("R",0)

SIR <= Reaction("infection",β, "S + I --> I + I")
SIR <= Reaction("recovery", γ, "I     --> R")

#SIRsim = simulate(SIR, Direct(), tfinal=days, rates_cache=HasRates, ntrials=50)
SIRsim = simulate(SIR, Direct(), tfinal=days, rates_cache=HasRates, save_points=0:1:365, ntrials=50)

plot(SIRsim, summary=:mean, label=["S" "I" "R"])

pSIR = plot(SIRsim[1], linecolor=["red" "blue" "green"],label=["S" "I" "R"])
for i=2:length(SIRsim)
    pSIR = plot!(SIRsim[i], summary = :trajectory, linecolor=["red" "blue" "green"], label="")
end
display(pSIR)

plot(SIRsim, timepoint=days, summary=:histogram, label=["S" "I" "R"],layout=(2,2),bins=20)

########################################################
### SIRD ################################################
########################################################
n_people = 1000
β₀       = 0.5
β        = β₀/n_people
γ        = 0.25
α        = 0.05
days     = 100

SIRD = Network("SIRD")
SIRD <= Species("S",n_people - 1)
SIRD <= Species("I",1)
SIRD <= Species("R",0)
SIRD <= Species("D",0)

SIRD <= Reaction("infection",β, "S + I --> I + I")
SIRD <= Reaction("recovery", γ, "I     --> R")
SIRD <= Reaction("death",α,     "I     --> D")

SIRDsim = simulate(SIRD, Direct(), tfinal=days, rates_cache=HasRates, ntrials=50)
#SIRsim = simulate(SIR, Direct(), tfinal=days, rates_cache=HasRates, save_points=0:1:365, ntrials=50)

plot(SIRDsim, summary=:mean, label=["S" "I" "R" "D"])

pSIRD = plot(SIRDsim[1], linecolor=["red" "blue" "green" "black"],label=["S" "I" "R" "D"])
for i=2:length(SIRDsim)
    pSIRD = plot!(SIRDsim[i], summary = :trajectory, linecolor=["red" "blue" "green" "black"], label="")
end

display(pSIRD)
