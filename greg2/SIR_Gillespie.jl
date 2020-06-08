using Gillespie

import Random: seed!

function SIRD_propensities(state, parameters)
    S, I, R, D = state
    N = S + I + R

    μ = parameters.μ
    β = parameters.β
    λ = parameters.λ
    α = parameters.α

    dS = μ * (1 - S) - beta * S * I / N
    dI = beta * S * I / N - (μ + λ + α) * I
    #dR = λ * I - μ * R
    dD = α * I

    return [dS, dI, dD]
end

transitions = [
               [ 1, -1,  0,  0]; # S -> I
               [ 0,  1, -1,  0]; # I -> R
               [-1,  0,  0,  1]; # S -> D
               [ 0, -1,  0,  1]; # I -> D
               [ 0,  0, -1,  1]; # R -> D
              ]

initial_state = [0.99, 0.01, 0.0]

parameters = (μ = 0.01,
              β = 0.3,
              λ = 0.05,
              α = 0.2)

end_time = 365

seed!(1234)

result = ssa(initial_state, SIRD_propensities, transitions, parameters, stop_time)

data = ssa_data(result)
