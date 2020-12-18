#function method
#%%
using DifferentialEquations
using Plots

# sir_model = @ode_def SIRModel begin
#     dS = -b*S*I
#     dI = b*S*I - g*I
#     dR = g*I
# end b g

p = (β = 0.3, γ = 0.1, σ=0.1)
u0 = [0.999, 0.0, 0.001, 0.0]
tspan = [0, 365.0] #1 year
tsave = 1.0:1.0:365
function sir_model(du, u, p, t)
    S, E, I, R = u
    b, g, s = p.β, p.γ, p.σ #b = transmission rate ; g = recovery rate
    du[1] = -b*S*I #note we are treating S,I,R as densities
    du[2] = b*S*I - g*E
    du[3] = g*E - s*I
    du[4] = s*I
end

problem = ODEProblem(sir_model, u0, tspan, p)
solution = solve(problem,saveat=tsave)
plot(solution, title = "SIR MODEL for COVID-19")
#%%
xx = solution[3,:]
plot(xx)

plot(log.(xx[2:365]./xx[1:364]))
