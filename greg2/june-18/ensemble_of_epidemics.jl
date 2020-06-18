include("stochastic_tools.jl")

n_people = 1000
β₀ = 0.5
γ₀ = 100

problem = stochastic_SIR_problem(n_people; β = β₀ / n_people, γ = γ₀ / n_people)

ensemble = solve_ensemble(problem, 100)

p = plot_solution(ensemble[1], alpha=0.1)

for i = 2:length(ensemble)
    plot_solution!(p, ensemble[i], alpha=0.1)
end

display(p)
