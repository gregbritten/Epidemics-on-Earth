using DifferentialEquations, Plots

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

exposure_rate(state, parameters, time) = transmission_rate(parameters.β, time) * state[1] * state[3]
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
                                time_span = (0, 100.0)
                                )

    # Parameters
    parameters = (β=β, σ=σ, γ=γ)

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

Returns t, S, I, R.
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

    r_plot = plot(t, R;
                  label = "R",
                  color = :black,
                  kwargs...)

    two_pane = plot(sei_plot, r_plot, layout=(2, 1))

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

    plot!(r_plot, t, R;
          label = false,
          color = :black,
          kwargs...)

    return nothing
end
