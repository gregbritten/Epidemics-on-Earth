using DifferentialEquations, Plots
using Pkg
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("OnlineStats")
Pkg.add("StatsPlots")
using StatsPlots
using DataFrames, CSV
using Statistics, OnlineStats
# SEIR problem...
#

# integrator.u[1] : "S"
# integrator.u[2] : "E"
# integrator.u[3] : "I"
# integrator.u[4] : "R"

# Define transitions for jump problem
function susceptible_to_exposed!(integrator)
    integrator.u[1] -= 1
    integrator.u[2] += 1
    return nothing
end

function exposed_to_infectious!(integrator)
    integrator.u[2] -= 1
    integrator.u[3] += 1
    return nothing
end

function infectious_to_resistant!(integrator)
    integrator.u[3] -= 1
    integrator.u[4] += 1
    return nothing
end

# Define rates for jump problem
#
# Below, we use multiple dispatch to obtain the transmission rate
# at the current time given the parameter "β".
# If `β` is a function, it is called with signature `β(time)`
# If β is anything else, it is simply returned.
transmission_rate(β, time) = β
transmission_rate(β::Function, time) = β(time)

exposure_rate(state, parameters, time) = transmission_rate(parameters.β, time) * state[1] * state[3] / parameters.n_people
infection_rate(state, parameters, time) = parameters.σ * state[2]
recovery_rate(state, parameters, time) = parameters.γ * state[3]

# To be strictly correct we need to use VariableRateJump here.
exposure = ConstantRateJump(exposure_rate, susceptible_to_exposed!)
infection = ConstantRateJump(infection_rate, exposed_to_infectious!)
recovery = ConstantRateJump(recovery_rate, infectious_to_resistant!)

"""
    stochastic_SIR_problem(n_people; β, γ, kwargs...)

Returns a `DiscreteProblem` that implements a stochastic SIR model
for the spread of a disease among `n_people` "well-mixed" individiuals
with recovery rate `γ` and `constant transmission rate `β` or time-varying
transmission rate `β(time)`.
"""
function stochastic_SIR_problem(n_people;
                                β,
                                σ,
                                γ,
                                percent_infectious = 1,
                                percent_exposed = 0,
                                time_span = (0, 300.0)
                                )

    # Parameters
    parameters = (β=β, σ=σ, γ=γ, n_people=n_people)

    # Initial condition
    n_infectious = ceil(n_people * percent_infectious / 100)
    n_exposed = ceil(n_people * percent_exposed / 100)
    n_susceptible = n_people - n_infectious - n_exposed

    # Make discrete problem
    discrete_problem = DiscreteProblem([n_susceptible, n_exposed, n_infectious, 0], time_span, parameters)

    return discrete_problem
end

"""
    solve_single(discrete_problem)

Wrap `discrete_problem` in a `JumpProblem` and solve it.
"""
function solve_single(discrete_problem)
    jump_problem = JumpProblem(discrete_problem, Direct(), exposure, infection, recovery)
    solution = solve(jump_problem, FunctionMap())
    print(solution)
    return solution
end

"""
    solve_ensemble(discrete_problem, n_ensemble=10)

Solve a `n_ensemble` `discrete_problem`s.
"""
function solve_ensemble(discrete_problem, n_ensemble=10)
    ensemble = []

    for i = 1:n_ensemble
        jump_problem = JumpProblem(discrete_problem, Direct(), exposure, infection, recovery)
        solution = solve(jump_problem, FunctionMap())
        push!(ensemble, solution)
    end

    return ensemble
end

"""
    unpack(solution)

Returns t, S, E, I, R.
"""
function unpack(solution)
    S = map(u -> u[1], solution.u)
    E = map(u -> u[2], solution.u)
    I = map(u -> u[3], solution.u)
    R = map(u -> u[4], solution.u)
    t = solution.t

    return t, S, E, I, R
end

function plot_solution(solution; kwargs...)
    t, S, E, I, R = unpack(solution)

    sei_plot = plot(t, S;
                    label = "S",
                    color = :blue,
                   xlabel = "Time (days)",
                   ylabel = "People",
                   kwargs...)

    plot!(sei_plot, t, E;
          label = "E",
          color = :red,
          kwargs...)

    plot!(sei_plot, t, I;
          label = "I",
          color = :black,
          kwargs...)
    # plot!(sei_plot, t, R;
    #     label = "R",
    #     color = :green,
    #     kwargs...)
    #
    r_plot = plot(t, R;
                  label = "R",
                  color = :black,
                  kwargs...)

    two_pane = plot(sei_plot, r_plot, layout=(2, 1))
    # two_pane = plot(sei_plot)

    return two_pane
end

function plot_solution!(two_pane, solution; kwargs...)
    t, S, E, I, R = unpack(solution)
    sei_plot = two_pane[1]
    r_plot = two_pane[2]

    plot!(sei_plot, t, S;
          label = false,
          color = :blue,
          kwargs...)

    plot!(sei_plot, t, E;
          label = false,
          color = :red,
          kwargs...)

    plot!(sei_plot, t, I;
          label = false,
          color = :black,
          kwargs...)

    # plot!(sei_plot, t, R;
    #     label = false,
    #     color = :green,
    #     kwargs...)

    plot!(r_plot, t, R;
          label = false,
          color = :black,
          kwargs...)

    return nothing
end



## Different Betas

#
# β(t) = 0.15 * (1 + cos(2π * t / 365)) #seasonal
# β(t) = .042*(3+cos(0.1π*t))
# plot(β)
#
# A(t) = (.75/sqrt(t))
# β(t) = A.(t)*(1+cos(0.2π*t))

## Running stochastic results and plotting

# sizes = [100,1000, 10000, 100000] #sizes for CSV files to test Systrom
# times = [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50]
# social_distancing_time = 15
# reopening_times = [45, 60, 75, 90, 105, 120 ,135, 150]
# reopening_time = 100
# β₀ = 0.4
# betas = [0.01, 0.2] #low and high beta for CSV files to test Systrom. Yields R_t of 1.1 and 2.0

## IGNORE HERE
#varying n statistics
# xvar1 = Vector() #n people
# yvar1 = Vector() #mean(R)/n people
# yvar2 = Vector() #std(R)/mean(R)
# yvar5 = Vector() #variance
# yvar6 = Vector() #skewness
# #varying social distancing statistics
# xvar2 = Vector() #social distancing time
# yvar3 = Vector() #mean(R)/n people
# yvar4 = Vector() #std(R)/mean(R)
# yvar7 = Vector() #variance
#
# yvar8 = Vector() #skewness
# #varying time between distancing and reopening statistics
# xvar3 = Vector() #time in between aka reopening - distancing times
# yvar9 = Vector() #proportion of resurgences
# yvar10 = Vector() #mean(R)/n people
# yvar11 = Vector() #std(R)/mean(R) = coefficient of variation
#
#
# step_β(time) = time < social_distancing_time ? β₀ / n_people : β₀ / 10 / n_people
#
# function distance_and_reopen(time) #full reopening
#     if time < social_distancing_time || time > reopening_time
#         return β₀/n_people
#     else
#         return β₀/10/n_people
#     end
# end
#
# function distance_and_reopen2(time) #partial reopening
#     if time < social_distancing_time
#         return β₀/n_people
#     elseif time > reopening_time
#         return β₀/5/n_people
#     else
#         return β₀/10/n_people
#     end
# end

# function reopen_when_low(time)
#     if time < social_distancing_time
#         return β₀/n_people
    # elseif

# step_β(time) = 0.15 * (1 + cos(2π * time / 365))
# plot(step_β)

# low = 0.04
# high = 0.1
# function step_β(time)
#     if 1 <= mod(time, 7) <= 4 #monday-thursday
#         return high
#     else #friday-sunday
#         return low
#     end
# end

#
# problem = stochastic_SIR_problem(n_people; β = distance_and_reopen, γ = 0.1, σ = 0.25)
# ensemble = solve_ensemble(problem, 75)

## CODE TO GENERATE 8 CSV FILES LOOK HERE
# sizes = [100, 1000, 10000, 100000]
sizes = [1000000]
for size in sizes
    #adjust beta and gamma as desired below
    problem = stochastic_SIR_problem(size, β = 1.0, γ = 0.1, σ = 0.1)
    ensemble = solve_ensemble(problem, 75) #75 realizations per ensemble
    alpha = 0.2
    p1 = plot_solution(ensemble[1], alpha=alpha) #plots single
    display(p1)

    reduced_ensemble = [] #option to add condition of how many are recovered
    for (i, member) in enumerate(ensemble)
        final_recovered = member.u[end][4]
        final_recovered > 0 && push!(reduced_ensemble, member)
    end

    p2 = plot_solution(reduced_ensemble[1], alpha=alpha)
    for i = 2:length(reduced_ensemble)
        plot_solution!(p2, reduced_ensemble[i], alpha=alpha) #what does this do
    end
    display(p2)

    ## Creating dataframe mapping days to number of new infections on that day
    fns = [daily_new_susceptible, daily_new_exposed, daily_new_infectious, daily_new_recovered]
    count = 1
    days = [1:300...]
    df = DataFrame(day = days)
    new_x = Vector{Float64}()
    for realization in ensemble
        for f in fns
            new_x = Vector{Float64}()
            for day in days
                new = f(realization, day)
                # print("NEWW", new)
                # print("LENGTh", length(new))
                push!(new_x, new)
            end
            title = count
            count += 1
            df[title] = new_x
        end

    end

    # df = DataFrame(day = days, new_inf = new_infections)
    CSV.write("C:\\Users\\spkho\\OneDrive\\Documents\\COVID UROP\\SEIR_n3_B3.csv",df)
end

# S1 = map(e -> e.u[1][1], ensemble) #susceptible at beginning
# S2 = map(e -> e.u[end][1], ensemble) #susceptible left
# R = map(e -> e.u[end][4], ensemble)
# print("S_start", S1)
# print("S_end", S2)
# print('R', R))

## General plotting for summer final pres
# display(p1)
# print(ensemble[2][2])
# R =  map(e -> e.u[end][4], ensemble) #number recovered = total infected
#
# reduced_ensemble = [] #option to add condition of how many are recovered
# for (i, member) in enumerate(ensemble)
#     final_recovered = member.u[end][4]
#     final_recovered > 0 && push!(reduced_ensemble, member)
# end
#
# p1 = plot_solution(reduced_ensemble[1], alpha=alpha)
# for i = 2:length(reduced_ensemble)
#     plot_solution!(p1, reduced_ensemble[i], alpha=alpha) #what does this do
# end
#
# plot!(p1, [1, 1] .* social_distancing_time, [0, n_people],
#       color = :green,
#       linewidth = 3,
#       alpha = alpha,
#       label = "Social distancing time")#adds label and line
#
# plot!(p1, [1, 1] .* reopening_time, [0, n_people],
#     color = :orange,
#     linewidth = 3,
#     alpha = alpha,
#     label = "Reopening time",
#     title = "Distancing and Partially Reopening")
# # display(p1)
#
# ##Probabilty density plot
# R = map(e -> e.u[end][4], ensemble) #number recovered = total infected
# print("mean = ", mean(R))
# print("variance = ", var(R))
# print("skewness = ", skewness(R))
# p2 = histogram(R / n_people, bin=100)
# display(p2)


function daily_new_susceptible(solution, day)
    """
    implement this if want to find decrease in susceptible, but that's
    just the opposite of daily new exposed
    """
    return 0
end

function daily_new_exposed(solution, day)
    """
    Calculating the number of new exposed persons on a specified day
    Compares change in susceptible
    """
    e_previous_day = searchsortedfirst(solution.t, day-1)
    previous_exposed = solution.u[e_previous_day][1]

    e_current_day = searchsortedfirst(solution.t, day)
    current_exposed = solution.u[e_current_day][1]

    change = previous_exposed - current_exposed

    if change >= 0
        return change
    else
        return 0.0
    end
end

function daily_new_infectious(solution, day)
    """
    calculating new number of infectious persons on a given day
    compares change in exposed
    """
    i_previous_day = searchsortedfirst(solution.t, day-1)
    previous_infectious = solution.u[i_previous_day][2]

    i_current_day = searchsortedfirst(solution.t, day)
    current_infectious = solution.u[i_current_day][2]

    change = previous_infectious - current_infectious

    if change >= 0
        return change
    else
        return 0.0
    end
end

function daily_new_recovered(solution, day)
    """
    calculating new number of recovered/"removed" persons on a given day
    compares change in infectious
    """
    r_previous_day = searchsortedfirst(solution.t, day-1)
    previous_recovered = solution.u[r_previous_day][3]

    r_current_day = searchsortedfirst(solution.t, day)
    current_recovered = solution.u[r_current_day][3]

    change = previous_recovered - current_recovered

    if change >= 0
        return change
    else
        return 0.0
    end
end


## Creating dataframe mapping days to number of new infections on that day
days = [1:100...]
new_infections = Vector{Float64}()
for day in days
    new = daily_new_infections(ensemble[2], day)
    push!(new_infections, new)
end
print(new_infections)

df = DataFrame(day = days, new_inf = new_infections)
# print(df)
## Export as CSV/Excel file
CSV.write("C:\\Users\\spkho\\OneDrive\\Documents\\COVID UROP\\test_large_population4.csv",df)

## Scatter plots: Finding trends in mean/std dev of total # infected based on population size / social distancing
#note to self: maybe use coefficient of variation instead of standard deviation to account for differnt pop size
# First: vary population size (n_people)
