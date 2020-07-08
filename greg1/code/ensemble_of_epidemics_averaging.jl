include("D:/Dropbox/Working/COVID19/UROP/github/greg2/june-18/stochastic_tools.jl")
using Statistics
using Interpolations
using StatsBase

n_people = 1000
β₀ = 0.5

problem = stochastic_SIR_problem(n_people; β = β₀/n_people, γ = 0.1)

ensemble = solve_ensemble(problem, 100)

p = plot_solution(ensemble[1], alpha=0.1)

for i = 2:length(ensemble)
    plot_solution!(p, ensemble[i], alpha=0.05)
end

display(p)

##--Extract and interpolate--#########
t, S, I, R = unpack(ensemble[2])

t_intp = 1.0:0.1:100
Sint   = LinearInterpolation(t,S)
Iint   = LinearInterpolation(t,I)
Rint   = LinearInterpolation(t,R)
sint   = Sint(t_intp)
iint   = Iint(t_intp)
rint   = Rint(t_intp)

plot(t_intp,[sint iint rint])

###################################

##--LOOP--############
E = Array{Float64}(undef,length(t_intp),3,100)
t_intp = 1.0:0.1:100

for i = 1:100
    t, S, I, R = unpack(ensemble[i])

    Sint   = LinearInterpolation(t,S)
    Iint   = LinearInterpolation(t,I)
    Rint   = LinearInterpolation(t,R)
    sint   = Sint(t_intp)
    iint   = Iint(t_intp)
    rint   = Rint(t_intp)

    E[:,:,i] = [sint iint rint]
end

E_mean = dropdims(mean(E,dims=3),dims=3)
E_sd   = dropdims(std(E,dims=3),dims=3)
E_med  = dropdims(median(E,dims=3),dims=3)
E_skew = dropdims(mapslices(skewness,E,dims=3),dims=3)
E_kurt = dropdims(mapslices(kurtosis,E,dims=3),dims=3)

plot(t_intp,E_mean)
plot(t_intp,E_med)
plot(t_intp,E_sd)
plot(t_intp,E_skew)
plot(t_intp,E_kurt)
