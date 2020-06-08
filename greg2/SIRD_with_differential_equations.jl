using DifferentialEquations, Plots

"""Compute the right-hand side of the SIRD equations."""
function compute_dSIRD_dt!(dSIRD_dt, SIRD, parameters, time)
    S, I, R, D = SIRD

    N = S + I + R

    μ = parameters.μ
    β = parameters.β(time)
    λ = parameters.λ
    α = parameters.α

    # Normalized by N₀:
    dSIRD_dt[1] = μ * (1 - S) - β * S * I / N
    dSIRD_dt[2] = β * S * I / N - (μ + λ + α) * I
    dSIRD_dt[3] = λ * I - μ * R
    dSIRD_dt[4] = α * I

    return nothing
end

# Set a problem with seasonally-varying transmission rate
year = 365.0
β(t) = 0.15 * (1 + cos(2π * t / year))

parameters = (μ = 0.002,
              β = β,
              λ = 0.05,
              α = 0.01)

timespan = (0, 4year)
initial_SIRD = [0.99, 0.01, 0.0, 0.0]

problem = ODEProblem(compute_dSIRD_dt!, initial_SIRD, timespan, parameters)

solution = solve(problem)

t = solution.t
S, I, R, D = [map(u -> u[i], solution.u) for i = 1:4]

states = plot(t, [S, I, R, D], label=["S" "I" "R" "D"],
              ylabel="Fraction of population")

transmission = plot(t, β.(t), label="β", xlabel="Days", ylabel="β")

plot(states, transmission, layout=(2, 1), size=(1000, 400))
