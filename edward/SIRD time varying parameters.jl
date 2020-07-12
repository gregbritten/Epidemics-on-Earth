using Plots

const niter = 7500
const dt = 0.1

const S = zeros(niter+1)
const I = zeros(niter+1)
const R = zeros(niter+1)
const D = zeros(niter+1)
const M = zeros(niter+1)

function sird_iter_model(μ, β, λ, α, N₀)
    S[1], I[1], R[1] = N₀

    for t in 1:niter
        dSdt   = -β(t)*S[t]*I[t]/(S[t]+I[t]+R[t]) + μ(t)*(1-S[t])
        dIdt   =  β(t)*S[t]*I[t]/(S[t]+I[t]+R[t]) - (μ(t)+λ(t)+α(t))*I[t]
        dRdt   =  λ(t)*I[t] - μ(t)*R[t]
        dDdt   =  α(t)*I[t]

        S[t+1] = S[t] + dSdt*dt
        I[t+1] = I[t] + dIdt*dt
        R[t+1] = R[t] + dRdt*dt
        D[t+1] = D[t] + dDdt*dt

        M[t+1] = dDdt*100
    end
end

β(t) = t < (0.2*niter) ? 0.02 : 0.06
λ(t) = 0.006
μ(t) = 0.01
α(t) = (I[t] > 0.25) ? (0.0025+0.05*(I[t]-0.25))/I[t] : 0.01

N = (0.99, 0.01, 0.0)
sird_iter_model(μ, β, λ, α, N)
plot((0:niter) .* dt, [S I R], label=['S' 'I' 'R'])
plot!((0:niter) .* dt, [M], label='M')

R_t(t) = β(t)/(λ(t)+μ(t)+α(t))
plot((0:niter) .* dt, [R_t(t) for t in (1:niter+1)])
