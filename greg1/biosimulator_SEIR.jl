
using BioSimulator
using Plots
using CSV
using DataFrames
########################################################
### SIRD ################################################
########################################################
n_people = 100000
β₀       = 0.3
β        = β₀/n_people
γ        = 0.1
σ        = 0.1
days     = 2*365

SEIR  = Network("SEIR")
SEIR <= Species("S",n_people - 1)
SEIR <= Species("E",0)
SEIR <= Species("I",1)
SEIR <= Species("R",0)

SEIR <= Reaction("exposure",  β, "S + I --> E + I")
SEIR <= Reaction("infection", σ, "E     --> I")
SEIR <= Reaction("recovery",  γ, "I     --> R")

tsave = 1.0:1.0:days

SEIRsim = simulate(SEIR, Direct(), tfinal=days, rates_cache=HasRates, ntrials=50,save_points=tsave)
#SIRsim = simulate(SIR, Direct(), tfinal=days, rates_cache=HasRates, save_points=0:1:365, ntrials=50)

#plot(SEIRsim, summary=:mean, label=["S" "E" "I" "R"])

pSEIR = plot(SEIRsim[1], linecolor=["red" "blue" "green" "black"],label=["S" "E" "I" "R"])
for i=2:length(SEIRsim)
    pSEIR = plot!(SEIRsim[i], summary = :trajectory, linecolor=["red" "blue" "green" "black"], label="")
end

display(pSEIR)

str = @sprintf "d:/dropbox/working/covid19/urop/github/greg1/SEIR_ensemble_n=%.0f_b=%.2f.csv" n_people β₀
CSV.write(str,SEIRsim)
