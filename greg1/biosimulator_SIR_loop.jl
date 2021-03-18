
using BioSimulator
using Plots
using CSV
using DataFrames
using Printf

########################################################
### SIRD ################################################
########################################################
#βs        = range(0.05,1.5,step=0.05)
#βs        = [range(0.01,0.04,step=0.01),range(0.05,1.5,step=0.05)]
βs         = 10 .^range(-2,0.2,step=0.05)
#n_peoples = [Int(1E3), Int(1E4), Int(1E5), Int(1E6)]
n_people = Int(1E5)
γ         = 0.1
σ         = 0.1     
R₀s       = βs/γ
days      = 365*5
ntrials   = Int(1E4)
tsave     = 1.0:1.0:days

for k=1:length(βs)
    β₀ = βs[k]
#    for j=1:4
#        n_people = n_peoples[j]
        β        = β₀/n_people

        SIR  = Network("SIR")
        SIR <= Species("S",n_people - 1)
        SIR <= Species("I",1)
        SIR <= Species("R",0)

        SIR <= Reaction("infection", β, "S + I --> I + I")
        SIR <= Reaction("recovery",  γ, "I     --> R")

function save_rates(simulator, state, model)
    copy(jump_rates(simulator))
end

    SIRsim = simulate(SIR, Direct(), tfinal=days, rates_cache=HasRates, ntrials=ntrials,
    save_points=tsave)
            #SEIRsim = simulate(SEIR, Direct(), tfinal=days, rates_cache=HasRates, ntrials=100,save_points=tsave,
            #    save_function =  save_rates)

        str = @sprintf "/users/gregorybritten/dropbox/working/covid19/urop/simulations_SIR/SIR_ensemble_n=%.0f_b=%.2f.csv" n_people β₀
        #str = @sprintf "d:/dropbox/working/covid19/urop/github/greg1/SEIR_ensemble_n=%.0f_b=%.2f_try.csv" n_people β₀
        CSV.write(str,SIRsim)
end
