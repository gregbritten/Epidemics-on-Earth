#function method
#%%
using DifferentialEquations
using Plots

# sir_model = @ode_def SIRModel begin
#     dS = -b*S*I
#     dI = b*S*I - g*I
#     dR = g*I
# end b g

p = (β = 0.1, γ = 0.1)
u0 = [0.99, 0.01, 0.0]
tspan = [0, 365.0] #1 year
function sir_model(du, u, p, t)
    S, I, R = u
    b, g = p.β, p.γ #b = transmission rate ; g = recovery rate
    du[1] = -b*S*I #note we are treating S,I,R as densities
    du[2] = b*S*I - g*I
    du[3] = g*I
end

problem = ODEProblem(sir_model, u0, tspan, p)
solution = solve(problem)
plot(solution, title = "SIR MODEL for COVID-19")
#%%
