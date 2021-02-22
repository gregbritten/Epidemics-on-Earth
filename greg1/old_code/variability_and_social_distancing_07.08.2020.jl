include("d:/Dropbox\\/Working\\/COVID19\\/UROP\\/github\\/greg1\\/stochastic_tools_07.08.2020.jl")
using Interpolations
using Statistics
using StatsBase

n_people = 10000
social_distancing_time = 30
β₀ = 0.5
step_β(time) = time < social_distancing_time ? β₀ / n_people : β₀ / 10 / n_people
#step_β(time) = 0.1*sin(time) + 1.0

problem = stochastic_SIR_problem(n_people; β = step_β, γ = 0.1, time_span = (0, 365.0), percent_infected = 1/n_people)

ensemble = solve_ensemble(problem, 200)

# Plot results
p = plot_solution(ensemble[1], alpha=0.2)


for i = 2:length(ensemble)
    plot_solution!(p, ensemble[i], alpha=0.2)
end

plot!(p, [1, 1] .* social_distancing_time, [0, n_people],
      color = :gray,
      linewidth = 3,
      alpha = 0.2,
      label = "Social distancing time")

display(p)


##--Extract and interpolate--#########
t, S, I, R = unpack(ensemble[2])

t_intp = 1.0:0.1:100.0
Sint   = LinearInterpolation(t,S)
Iint   = LinearInterpolation(t,I)
Rint   = LinearInterpolation(t,R)
sint   = Sint(t_intp)
iint   = Iint(t_intp)
rint   = Rint(t_intp)

plot(t_intp,[sint iint rint])

###################################

##--LOOP--############
E = Array{Float64}(undef,length(t_intp),3,200)
t_intp = 1.0:0.1:100

for i = 1:200
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

p1 = plot(t_intp,E_mean,label=["S" "I" "R"],title="Mean")
#p2 = plot(t_intp,E_med,label=["S" "I" "R"])
p3 = plot(t_intp,E_sd,label=["S" "I" "R"],title="Standard Deviation")
p4 = plot(t_intp,E_sd./E_mean,label=["S" "I" "R"],title="SD/Mean")
p5 = plot(t_intp,E_skew,label=["S" "I" "R"],title="Skewness")
plot(p1,p3,p4,p5)

Ecs = cumsum(E,dims=1)

hS  = histogram(E[length(t_intp),1,:],bins=30,title="Susceptible",label="")
hI  = histogram(E[length(t_intp),2,:],bins=30,title="Infected",label="")
hIc = histogram(Ecs[length(t_intp),2,:],bins=30,title="Cumulative Infected",label="")
hR  = histogram(E[length(t_intp),3,:],bins=30,title="Removed",label="")
plot(hS,hI,hIc,hR)
