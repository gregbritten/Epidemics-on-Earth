using DifferentialEquations, Plots

# This script sets up a stochastic SIR model.
#
# See https://docs.sciml.ai/stable/tutorials/discrete_stochastic_example/
#
# for inspiration.

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

# Define rates
function infection_rate(state, parameters, time)
    S = state[1]
    I = state[2]
    return parameters.β * S * I
end

function recovery_rate(state, parameters, time)
    I = state[2]
    return parameters.γ * I
end

infection = ConstantRateJump(infection_rate, susceptible_to_infectious!)
recovery = ConstantRateJump(recovery_rate, infectious_to_resistant!)

# Parameters
n_people = 100

parameters = (
    β = 0.5 / n_people,
    γ = 0.1,
)

time_span = (0, 100.0)

# Initial condition
n_infected = 1
n_susceptible = n_people - n_infected

# Make discrete problem
discrete_problem = DiscreteProblem([n_susceptible, n_infected, 0], time_span, parameters)

jump_problem = JumpProblem(discrete_problem, Direct(), infection, recovery)
solution = solve(jump_problem, FunctionMap())

S = map(u -> u[1], solution.u)
I = map(u -> u[2], solution.u)
R = map(u -> u[3], solution.u)
t = solution.t

si_plot = plot(t, S;
                label = "S",
                color = :blue,
               xlabel = "Time (days)",
               ylabel = "People")

plot!(si_plot, t, I;
      label = "I",
      color = :red)

r_plot = plot(t, R;
              label = "R",
              color = :black)

two_pane = plot(si_plot, r_plot, layout=(2, 1))

display(two_pane)
