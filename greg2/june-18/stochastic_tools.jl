using DifferentialEquations, Plots

# Define transitions for jump problem
function susceptible_to_infectious!(integrator)
    integrator.u[1] -= 1
    integrator.u[2] += 1
    return nothing
end

function infectious_to_resistant!(integrator)
    integrator.u[2] -= 1
    integrator.u[3] += 1
    return nothing
end

# Define rates for jump problem
transmission_rate(β::Number, time) = β
transmission_rate(β::Function, time) = β(time)

infection_rate(state, parameters, time) = transmission_rate(parameters.β, time) * state[1] * state[2]
recovery_rate(state, parameters, time) = parameters.γ * state[2]

infection = ConstantRateJump(infection_rate, susceptible_to_infectious!)
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
                                γ,
                                percent_infected = 1,
                                time_span = (0, 100.0)
                                )

    # Parameters
    parameters = (β=β, γ=γ)

    # Initial condition
    n_infected = ceil(n_people * percent_infected / 100)
    n_susceptible = n_people - n_infected

    # Make discrete problem
    discrete_problem = DiscreteProblem([n_susceptible, n_infected, 0], time_span, parameters)

    return discrete_problem
end

function solve_single(discrete_problem)
    jump_problem = JumpProblem(discrete_problem, Direct(), infection, recovery)
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
        jump_problem = JumpProblem(discrete_problem, Direct(), infection, recovery)
        solution = solve(jump_problem, FunctionMap())
        push!(ensemble, solution)
    end

    return ensemble
end

function unpack(solution)
    S = map(u -> u[1], solution.u)
    I = map(u -> u[2], solution.u)
    R = map(u -> u[3], solution.u)
    t = solution.t

    return t, S, I, R
end

function plot_solution(solution; kwargs...)
    t, S, I, R = unpack(solution)

    si_plot = plot(t, S;
                    label = "S",
                    color = :blue,
                   xlabel = "Time (days)",
                   ylabel = "People",
                   kwargs...)

    plot!(si_plot, t, I;
          label = "I",
          color = :red,
          kwargs...)

    r_plot = plot(t, R;
                  label = "R",
                  color = :black,
                  kwargs...)

    two_pane = plot(si_plot, r_plot, layout=(2, 1))

    return two_pane
end

function plot_solution!(two_pane, solution; kwargs...)
    t, S, I, R = unpack(solution)

    si_plot = two_pane[1]
    r_plot = two_pane[2]

    plot!(si_plot, t, S;
          label = false,
          color = :blue,
          kwargs...)

    plot!(si_plot, t, I;
          label = false,
          color = :red,
          kwargs...)

    plot!(r_plot, t, R;
          label = false,
          color = :black,
          kwargs...)

    return nothing
end
