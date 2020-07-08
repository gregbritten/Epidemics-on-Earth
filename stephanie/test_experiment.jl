include("SEIR.jl")

n_people = 100
social_distancing_time = 15
β₀ = 0.5
step_β(time) = time < social_distancing_time ? β₀ / n_people : β₀ / 10 / n_people

problem = stochastic_SIR_problem(n_people; β = step_β, γ = 0.1, σ = 0.25)

ensemble = solve_ensemble(problem, 100)

# Plot results
alpha = 0.2
p1 = plot_solution(ensemble[1], alpha=alpha)

for i = 2:length(ensemble)
    plot_solution!(p1, ensemble[i], alpha=alpha)
end

plot!(p1, [1, 1] .* social_distancing_time, [0, n_people],
      color = :gray,
      linewidth = 3,
      alpha = alpha,
      label = "Social distancing time")

display(p1)

R = map(e -> e.u[end][4], ensemble)

p2 = histogram(R / n_people, bin=10)

display(p2)
