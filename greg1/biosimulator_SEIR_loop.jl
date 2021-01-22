
using BioSimulator
using Plots
using CSV
using DataFrames
using Printf

########################################################
### SIRD ################################################
########################################################
βs        = [0.1, 0.15, 0.2, 0.25, 0.3]
n_peoples = [Int(1E3), Int(1E4), Int(1E5), Int(1E6)]
γ         = 0.1
σ         = 0.1
days      = 365*5
ntrials   = Int(1E3)
tsave     = 1.0:1.0:days

for k=1:5
    β₀ = βs[k]
    for j=1:4
        n_people = n_peoples[j]
        β        = β₀/n_people

        SEIR  = Network("SEIR")
        SEIR <= Species("S",n_people - 1)
        SEIR <= Species("E",0)
        SEIR <= Species("I",1)
        SEIR <= Species("R",0)

        SEIR <= Reaction("exposure",  β, "S + I --> E + I")
        SEIR <= Reaction("infection", σ, "E     --> I")
        SEIR <= Reaction("recovery",  γ, "I     --> R")

function save_rates(simulator, state, model)
    copy(jump_rates(simulator))
end

        #SEIRsim = simulate(SEIR, Direct(), tfinal=days, rates_cache=HasRates, ntrials=ntrials,save_points=tsave,
    #        save_function = save_rates)
            SEIRsim = simulate(SEIR, Direct(), tfinal=days, rates_cache=HasRates, ntrials=100,save_points=tsave, 
                save_function =  save_rates)

        str = @sprintf "d:/dropbox/working/covid19/urop/github/greg1/SEIR_ensemble_n=%.0f_b=%.2f.csv" n_people β₀
        str = @sprintf "d:/dropbox/working/covid19/urop/github/greg1/SEIR_ensemble_n=%.0f_b=%.2f_try.csv" n_people β₀
        CSV.write(str,SEIRsim)
    end
end
