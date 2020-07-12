using DifferentialEquations
using Plots

const niter = 8000
const dt = 0.1

const S = zeros(niter+1)
const I = zeros(niter+1)
const R = zeros(niter+1)
const D = zeros(niter+1)

function sird_iter_model(μ, β, λ, α, N₀)
    S[1], I[1], R[1], _ = N₀

    for t in 1:niter
        dSdt   = -β*S[t]*I[t]/(S[t]+I[t]+R[t]) + μ*(1-S[t])
        dIdt   =  β*S[t]*I[t]/(S[t]+I[t]+R[t]) - (μ+λ+α)*I[t]
        dRdt   =  λ*I[t] - μ*R[t]
        dDdt   =  α*I[t]

        S[t+1] = S[t] + dSdt*dt
        I[t+1] = I[t] + dIdt*dt
        R[t+1] = R[t] + dRdt*dt
        D[t+1] = D[t] + dDdt*dt
    end
end

function sird_dfeq_model(μ, β, λ, α, N₀)
    function sird_ode(du, u, p, t)
        S, I, R, D = u
        μ, β, λ, α = p
        du[1] = -β*S*I/(S+I+R) + μ*(1-S)
        du[2] =  β*S*I/(S+I+R) - (μ+λ+α)*I
        du[3] = λ*I - μ*R
        du[4] = α*I
    end

    sird_prob_func = ODEProblem(sird_ode, N₀, (0, niter*dt), (μ, β, λ, α))
    solve(sird_prob_func)
end

function sird_steady_test(μ, β, λ, α)
    R₀ = β/(λ+μ+α)
    if R₀ > 1.0
        N∞ = (β-α*R₀)/(β-α)
        N∞/R₀, μ*(R₀-N∞)/β, λ*(R₀-N∞)/β
    else
        1.0, 0.0, 0.0
    end
end

# Epidemic model
μ, β, λ, α = 0.01, 0.2, 0.006, 0.02
N₀ = [0.99, 0.01, 0.0, 0.0]

# Healthy equilibrium R₀ < 1
μ, β, λ, α = 0.01, 0.006, 0.2, 0.02
N₀ = [0.33, 0.33, 0.34, 0.0]

# R₀ = 1
μ, β, λ, α = 0.01, (0.01+0.006+0.02), 0.006, 0.02
N₀ = [0.9, 0.1, 0.0, 0.0]

# Iterative
sird_iter_model(μ, β, λ, α, N₀)
plot((0:niter) .* dt, [S I R D], label=['S' 'I' 'R' 'D'])
# ODE
ode_sol = sird_dfeq_model(μ, β, λ, α, N₀)
plot(ode_sol)

# Prediction
S∞, I∞, R∞ = sird_steady_test(μ, β, λ, α)
hline!([S∞, I∞, R∞])
