using DifferentialEquations, Plots

# SEIR problem...
#
# integrator.u[1] : "S"
# integrator.u[2] : "E"
# integrator.u[3] : "M"
# integrator.u[4] : "V"
# integrator.u[5] : "I"
# integrator.u[6] : "R"

# Define transitions for jump problem
function susceptible_to_exposed!(integrator)
    integrator.u[1] -= 1
    integrator.u[2] += 1
    return nothing
end

function exposed_to_medium_infectious!(integrator)
    integrator.u[2] -= 1
    integrator.u[3] += 1
    return nothing
end

function medium_infectious_to_very_infectious!(integrator)
    integrator.u[3] -= 1
    integrator.u[4] += 1
    return nothing
end

function very_infectious_to_slightly_infectious!(integrator)
    integrator.u[4] -= 1
    integrator.u[5] += 1
    return nothing
end

function slightly_infectious_to_resistant!(integrator)
    integrator.u[5] -= 1
    integrator.u[6] += 1
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

exposure_rate(state, parameters, time) = (
                                            transmission_rate(parameters.βᴹ, time) * state[1] * state[3]
                                          + transmission_rate(parameters.βⱽ, time) * state[1] * state[4]
                                          + transmission_rate(parameters.βᴵ, time) * state[1] * state[5]
                                         )

medium_infectious_rate(state, parameters, time)   = parameters.σᴹ * state[2]
very_infectious_rate(state, parameters, time)     = parameters.σⱽ * state[3]
slightly_infectious_rate(state, parameters, time) = parameters.σᴵ * state[4]

recovery_rate(state, parameters, time) = parameters.γ * state[5]

# To be strictly correct we need to use VariableRateJump here.
exposure = ConstantRateJump(exposure_rate, susceptible_to_exposed!)

medium_infection   = ConstantRateJump(medium_infectious_rate,   exposed_to_medium_infectious!)
very_infection     = ConstantRateJump(very_infectious_rate,     medium_infectious_to_very_infectious!)
slightly_infection = ConstantRateJump(slightly_infectious_rate, very_infectious_to_slightly_infectious!)

recovery = ConstantRateJump(recovery_rate, slightly_infectious_to_resistant!)

"""
    stochastic_SIR_problem(n_people; β, γ, kwargs...)

Returns a `DiscreteProblem` that implements a stochastic SIR model
for the spread of a disease among `n_people` "well-mixed" individiuals
with recovery rate `γ` and `constant transmission rate `β` or time-varying
transmission rate `β(time)`.
"""
function stochastic_SIR_problem(n_people;
                                βᴹ,
                                βⱽ,
                                βᴵ,
                                σᴹ,
                                σⱽ,
                                σᴵ,
                                γ,
                                percent_infectious = 1,
                                percent_exposed = 0,
                                time_span = (0, 100.0)
                                )

    # Parameters
    parameters = (
                  βᴹ = βᴹ, 
                  βⱽ = βⱽ, 
                  βᴵ = βᴵ, 
                  σᴹ = σᴹ, 
                  σⱽ = σⱽ, 
                  σᴵ = σᴵ, 
                  γ = γ
                 )

    # Initial condition
    n_infectious = ceil(n_people * percent_infectious / 100)
    n_exposed = ceil(n_people * percent_exposed / 100)
    n_susceptible = n_people - n_infectious - n_exposed

    # Make discrete problem
    discrete_problem = DiscreteProblem([n_susceptible, n_exposed, n_infectious, 0, 0, 0], time_span, parameters)

    return discrete_problem
end

"""
    solve_single(discrete_problem)

Wrap `discrete_problem` in a `JumpProblem` and solve it.
"""
function solve_single(discrete_problem)
    jump_problem = JumpProblem(discrete_problem, Direct(), exposure, medium_infection, very_infection, slightly_infection, recovery)
    solution = solve(jump_problem, FunctionMap())
    return solution
end

"""
    solve_ensemble(discrete_problem, n_ensemble=10)

Solve a `n_ensemble` `discrete_problem`s.
"""
function solve_ensemble(discrete_problem, n_ensemble=10)
    ensemble = []

    for i = 1:n_ensemble
        jump_problem = JumpProblem(discrete_problem, Direct(), exposure, medium_infection, very_infection, slightly_infection, recovery)
        solution = solve(jump_problem, FunctionMap())
        push!(ensemble, solution)
    end

    return ensemble
end

"""
    unpack(solution)

Returns t, S, I, R.
"""
function unpack(solution)
    S = map(u -> u[1], solution.u)
    E = map(u -> u[2], solution.u)
    M = map(u -> u[3], solution.u)
    V = map(u -> u[4], solution.u)
    I = map(u -> u[5], solution.u)
    R = map(u -> u[6], solution.u)
    t = solution.t

    return t, S, E, M, V, I, R
end

function plot_solution(solution; kwargs...)
    t, S, E, M, V, I, R = unpack(solution)

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

    plot!(sei_plot, t, M;
          label = "M",
          color = :black,
          kwargs...)

    plot!(sei_plot, t, V;
          label = "V",
          color = :purple,
          kwargs...)

    plot!(sei_plot, t, I;
          label = "I",
          color = :green,
          kwargs...)

    r_plot = plot(t, R;
                  label = "R",
                  color = :black,
                  kwargs...)

    two_pane = plot(sei_plot, r_plot, layout=(2, 1))

    return two_pane
end

function plot_solution!(two_pane, solution; kwargs...)
    t, S, E, M, V, I, R = unpack(solution)

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

    plot!(sei_plot, t, M;
          label = false,
          color = :black,
          kwargs...)

    plot!(sei_plot, t, V;
          label = false,
          color = :purple,
          kwargs...)

    plot!(sei_plot, t, I;
          label = false,
          color = :green,
          kwargs...)

    plot!(r_plot, t, R;
          label = false,
          color = :black,
          kwargs...)

    return nothing
end
