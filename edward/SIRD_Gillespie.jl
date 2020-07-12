using Random
using Plots

niter = Int(3E3)

# Reactions
# 1. S + I -> 2I
# 2. I -> R
# 3. I -> D

const Ss = zeros(niter+1)
const Is = zeros(niter+1)
const Rs = zeros(niter+1)
const Ds = zeros(niter+1)

const Ts = zeros(niter+1)

# Initial conditions
N = 1000
I₀ = 1
const x = [N-I₀, I₀, 0, 0]
Ss[1], Is[1] = x[1], x[2]

# Parameters
μ = 0.01
β, λ, α = 0.07, 0.006, 0.02

# infection, death, recovery
a_i(x) = β*x[1]*x[2]/(x[1]+x[2]+x[3])
a_d(x) = α*x[2]
a_r(x) = λ*x[2]
a_0(x) = a_i(x) + a_d(x) + a_r(x)

for i in 1:niter
    r₁ = rand(Float64)
    r₂ = rand(Float64)
    τ  = -log(r₁)/a_0(x)

    # infection
    if a_i(x) > r₂*a_0(x)
        x[1] -= 1
        x[2] += 1
    # death
    elseif a_i(x)+a_d(x) > r₂*a_0(x)
        x[2] -= 1
        x[4] += 1
    # recovery
    elseif r₂ < 1
        x[2] -= 1
        x[3] += 1
    end

    Ss[i+1] = x[1]
    Is[i+1] = x[2]
    Rs[i+1] = x[3]
    Ds[i+1] = x[4]
    Ts[i+1] = τ
end

t = cumsum(Ts)
plot(t, [Ss Is Rs Ds],linetype=:steppost)
