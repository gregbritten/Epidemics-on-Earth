include("stochastic_tools.jl")

n_people = 1000
social_distancing_time = 10
β₀ = 0.5
step_β(time) = time < social_distancing_time ? β₀ / n_people : β₀ / 10 / n_people

problem = stochastic_SIR_problem(n_people; β = step_β, γ = 0.1)

#problem = stochastic_SIR_problem(n_people; β = β₀ / n_people, γ = 0.1)

ensemble = solve_ensemble(problem, 10)

# Plot results
alpha = 0.2
p = plot_solution(ensemble[1], alpha=alpha)

for i = 2:length(ensemble)
    plot_solution!(p, ensemble[i], alpha=alpha)
end

plot!(p, [1, 1] .* social_distancing_time, [0, n_people],
      color = :gray,
      linewidth = 3,
      alpha = alpha,
      label = "Social distancing time")

display(p)
