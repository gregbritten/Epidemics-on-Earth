using DifferentialEquations
using Plots

dt = 0.1
const t_vacc = 85            #time until vaccination is introduced in days
T = .5                       #period between periodic vaccinations
global niter = Int((365)/dt)
global n_1iter = Int64((t_vacc)/dt)
global Titer = Int((365-t_vacc)/T)


S = zeros(niter+1)
E = zeros(niter+1)
I = zeros(niter+1)
R = zeros(niter+1)
D = zeros(niter+1)
I[1] = .01
S[1] = 1 - S[1]

# parameters in the PVS SEIR process
β = .1                #Transmission Rate
α = .1               #Rate E -> I
γ = .05               #Recovery Rate
λ = .001              #Mortality Rate
p = .1                #Proportion of Suspectibles Vaccinated in each Pulse


for t in 1:(n_1iter - 1)
    dSdt = -λ*S[t] - β*S[t]*I[t]
    dEdt = β*S[t]*I[t] - (α+λ)*E[t]
    dIdt = α*E[t] - (γ+λ)*I[t]
    dRdt = γ*I[t] - λ*R[t]
    dDdt = λ*(S[t]+I[t]+E[t]+R[t])

    S[t+1] = S[t] + dSdt*dt
    E[t+1] = E[t] + dEdt*dt
    I[t+1] = I[t] + dIdt*dt
    R[t+1] = R[t] + dRdt*dt
    D[t+1] = D[t] + dDdt*dt
end

for t in n_1iter:niter
    if (t-t_vacc)/T == 0
        for i in 0:Titer
            S[t_vacc+i*T+1] = S[t_vacc + i*T] *(1-p)
        end
        continue
    else
        while((t-t_vacc)/T != 0)
            S[t+1] = S[t] + dSdt*dt
            E[t+1] = E[t] + dEdt*dt
            I[t+1] = I[t] + dIdt*dt
            R[t+1] = R[t] + dRdt*dt
            D[t+1] = D[t] + dDdt*dt
        end
        continue
    end
end
end
println("")
#plot([S E I R D], label=['S' 'E' 'I' 'R' 'D'])
