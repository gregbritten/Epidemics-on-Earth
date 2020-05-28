using DifferentialEquations
using Plots

sir_ode = @ode_def SIRModel begin
    dS = -β*S*I
    dI = β*S*I-γ*I
    dR = γ*I
end β γ

parms = [0.1,0.05]
init  = [0.99,0.01,0.0]
tspan = (0.0,200.0)

sir_prob = ODEProblem(sir_ode,init,tspan,parms)

sir_sol = solve(sir_prob,saveat = 0.1);

plot(sir_sol,xlabel="Time",ylabel="Number")

##################################################
# AS A FUNCTION ##################################

function sir_ode(du,u,p,t)
    S,I,R = u
    b,g   = p
    du[1] = -b*S*I
    du[2] = b*S*I-g*I
    du[3] = g*I
end

sir_prob_func = ODEProblem(sir_ode2,init,tspan,parms)
sir_sol_func  = solve(sir_prob2,saveat = 0.1)

plot(sir_sol_func)

##################################################
dt    = 0.1
niter = Int(365/dt)

S = zeros(niter+1)
I = zeros(niter+1)
R = zeros(niter+1)
I[1] = 0.01
S[1] = 1 - I[1]

β = 0.1
γ = 0.05

for t in 1:niter
    dSdt   = -β*S[t]*I[t]
    dIdt   =  β*S[t]*I[t] - γ*I[t]
    dRdt   =  γ*I[t]

    S[t+1] = S[t] + dSdt*dt
    I[t+1] = I[t] + dIdt*dt
    R[t+1] = R[t] + dRdt*dt
end

plot([S I R], label=['S' 'I' 'R'])
