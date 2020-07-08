include("SEIR.jl")

function simulate_ensemble(n_people, social_distancing_time=15, β₀=0.5)

    step_β(time) = time < social_distancing_time ? β₀ / n_people : β₀ / 10 / n_people

    problem = stochastic_SIR_problem(n_people; β = step_β, γ = 0.1, σ = 0.25)

    ensemble = solve_ensemble(problem, 1000)

    R = map(e -> e.u[end][4], ensemble)

    p = histogram(R / n_people, bin=40, xlim=(0, 1))

    return p
end

p1 = simulate_ensemble(100)
p2 = simulate_ensemble(200)
p3 = simulate_ensemble(400)

plot(p1, p2, p3, layout=(3, 1))
