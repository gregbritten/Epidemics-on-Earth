using Plots

##################################################
dt    = 0.1
niter = Int(365/dt)
T     = range(0, stop=365 , length=Int(365/dt)+1)

S = zeros(niter+1)
I = zeros(niter+1)
R = zeros(niter+1)
D = zeros(niter+1)

I[1] = 0.01
S[1] = 1.0 - I[1]

μ = 0.01 #birth/death [/time]
β = 0.1 #tranmission rate [/time]
λ = 0.05 #recovery rate [/time]
α = 0.01 #death rate from infection [/time]
γ = 0.01 #vaccination rate [/time]

for t in 1:niter
    dSdt = μ - μ*S[t] - β*S[t]*I[t]/(S[t] + I[t] + R[t]) - γ*S[t]
    dIdt = β*S[t]*I[t]/(S[t] + I[t] + R[t]) - (μ + λ + α)*I[t]
    dRdt = λ*I[t] - μ*R[t] + γ*S[t]
    dDdt = α*I[t]

    S[t+1] = S[t] + dSdt*dt
    I[t+1] = I[t] + dIdt*dt
    R[t+1] = R[t] + dRdt*dt
    D[t+1] = D[t] + dDdt*dt
end

plot(T,[S I R D], label=['S' 'I' 'R' 'D'],
    xlabel="Days", ylabel="Population Fraction")
