include("SEIR.jl")

n_people = 100
social_distancing_time = 15
β₀ = 0.5

step_β(time) = time < social_distancing_time ? β₀ / n_people : β₀ / 10 / n_people

problem = stochastic_SIR_problem(n_people; β = step_β, γ = 0.1, σ = 0.25)

solution = solve_single(problem)

# Plot results
alpha = 0.8

p1 = plot_solution(solution, alpha=alpha)

plot_solution!(p1, solution, alpha=alpha)

plot!(p1, [1, 1] .* social_distancing_time, [0, n_people],
      color = :gray,
      linewidth = 3,
      alpha = alpha,
      label = "Social distancing time")

display(p1)

function new_infections(solution, day)
    i_previous_day = searchsortedfirst(solution.t, day-1)
    previous_infections = solution.u[i_previous_day][3]

    i_current_day = searchsortedfirst(solution.t, day)
    current_infections = solution.u[i_current_day][3]

    return current_infections - previous_infections
end

