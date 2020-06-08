using Plots, DifferentialEquations: ODEProblem

"""Compute the right-hand side of the SIRD equations."""
function compute_dSIRD_dt!(dSIRD_dt, SIRD, parameters, time)
    dS, dI, dR, dD = dSIRD_dt
    S,   I,  R,  D = SIRD

    N = S + I + R

    μ = parameters.μ
    β = parameters.β(time)
    λ = parameters.λ
    α = parameters.α

    # Normalized by N₀:
    dSdt = μ * (1 - S) - β * S * I / N
    dIdt = β * S * I / N - (μ + λ + α) * I
    dRdt = λ * I - μ * R
    dDdt = α * I

    return nothing
end

# Set a problem with time-dependent transmission rate:
year = 365
β₀ = 6e-6
β(t) = β₀ * (t - year / 2)^2 + 0.1

parameters = (μ = 0.01,
              β = β,
              λ = 0.05,
              α = 0.01)

timespan = (0, year)

problem = ODEProblem(compute_dSIRD_dt!, timespan, parameters)

solution = solve(problem)

t = solution.t
S, I, R, D = [map(u -> u[i], solution.u) for i = 1:4]

plot(t, [S I R D], label=["S", "I", "R", "D"], xlabel="Days", ylabel="Fraction of population")
