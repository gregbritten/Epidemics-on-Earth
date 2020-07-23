using Random, Plots

struct My_Parameters
    μ; β; a; γ; α
end

const N_SIM = Int(1e4)
const N_COMPARTMENT = 5

const u = zeros(Int, N_COMPARTMENT, N_SIM+1)
const T = zeros(N_SIM+1)

function birth!(u, i)
    u[1, i] += 1
end

function susceptible_to_exposed!(u, i)
    u[1, i] -= 1
    u[2, i] += 1
end

function exposed_to_infectious!(u, i)
    u[2, i] -= 1
    u[3, i] += 1
end

function infectious_to_resistant!(u, i)
    u[3, i] -= 1
    u[4, i] += 1
end

# Death events
function susceptible_death!(u, i)
    u[1, i] -= 1
    return nothing
end

function exposed_death!(u, i)
    u[2, i] -= 1
end

function infected_death!(u, i)
    u[3, i] -= 1
    u[5, i] += 1
end

function recovered_death!(u, i)
    u[4, i] -= 1
end

# Rates
birth_rate(state, parameters) = parameters.μ * sum(state)

exposure_rate(state, parameters) = parameters.β * state[1] * state[3] / (sum(state[1:4]))
infection_rate(state, parameters) = parameters.a * state[2]
recovery_rate(state, parameters) = parameters.γ * state[3]

susceptible_death_rate(state, parameters) = parameters.μ * state[1]
exposed_death_rate(state, parameters) = parameters.μ * state[2]
infected_death_rate(state, parameters) = (parameters.α+parameters.μ) * state[3]
recovered_death_rate(state, parameters) = parameters.μ * state[4]


function ∑a(state, parameters, i)
    a_i =  [birth_rate(state, parameters),

            exposure_rate(state, parameters),
            infection_rate(state, parameters),
            recovery_rate(state, parameters),

            susceptible_death_rate(state, parameters),
            exposed_death_rate(state, parameters),
            infected_death_rate(state, parameters),
            recovered_death_rate(state, parameters),
            ]
    # sort!(a_i)
    # @show a_i
    return sum(a_i[1:i])
end

function direct_ssa(initial_condition, parameters)
    u[:, 1] .= initial_condition

    for i in 1:N_SIM
        state = u[:, i]
        u[:, i+1] .= state

        a₀ = ∑a(state, parameters, 8)

        r₁ = rand(Float64)
        r₂ = rand(Float64)
        τ  = -log(r₁)/a₀

        # birth
        if r₂*a₀ <= ∑a(state, parameters, 1)
            birth!(u, i+1)
        # exposure
        elseif r₂*a₀ <= ∑a(state, parameters, 2)
            susceptible_to_exposed!(u, i+1)
        # infection
        elseif r₂*a₀ <= ∑a(state, parameters, 3)
            exposed_to_infectious!(u, i+1)
        # recovery
        elseif r₂*a₀ <= ∑a(state, parameters, 4)
            infectious_to_resistant!(u, i+1)
        elseif r₂*a₀ <= ∑a(state, parameters, 5)
            susceptible_death!(u, i+1)
        elseif r₂*a₀ <= ∑a(state, parameters, 6)
            exposed_death!(u, i+1)
        elseif r₂*a₀ <= ∑a(state, parameters, 7)
            infected_death!(u, i+1)
        else
            recovered_death!(u, i+1)
        end

        T[i+1] = τ
    end
end

function unpack(u)
    S = u[1, :]
    E = u[2, :]
    I = u[3, :]
    R = u[4, :]
    D = u[5, :]
    t = cumsum(T)
    return t, S, E, I, R, D
end

function plot_solution(u)
    t, S, E, I, R, D = unpack(u)
    return plot(t, [S E I R],
                label=['S' 'E' 'I' 'R'],
                xlabel = "Time",
                ylabel = "People",
                title="SEIRD Stochastic Simulation",
                linetype=:steppost)
end

μ=0.01; β=0.07; a=0.07; γ=0.006; α=0.02;
parameters = My_Parameters(μ, β, a, γ, α)
initial_condition = 950, 0, 50, 0, 0
direct_ssa(initial_condition, parameters)
plot_solution(u)
