using Random, Plots

struct My_Parameters
    μ; β; a; γ; α
end

const N_SIM = Int(3e3)
const N_COMPARTMENT = 5

const u = zeros(Int, N_COMPARTMENT, N_SIM+1)
const T = zeros(N_SIM+1)

function birth!(u, i)
    if (u[1, i] > 0) u[1, i] -= 1 end
    return nothing
end

function susceptible_to_exposed!(u, i)
    if (u[1, i] > 0)
        u[1, i] -= 1
        u[2, i] += 1
    end
    return nothing
end

function exposed_to_infectious!(u, i)
    if (u[2, i] > 0)
        u[2, i] -= 1
        u[3, i] += 1
    end
    return nothing
end

function infectious_to_resistant!(u, i)
    if (u[3, i] > 0)
        u[3, i] -= 1
        u[4, i] += 1
    end
    return nothing
end

# Death events
function exposed_death!(u, i)
    if (u[2, i] > 0)
        u[2, i] -= 1
        u[5, i] += 1
    end
    return nothing
end

function infected_death!(u, i)
    if (u[3, i] > 0)
        u[3, i] -= 1
        u[5, i] += 1
    end
    return nothing
end

function recovered_death!(u, i)
    if (u[4, i] > 0)
        u[4, i] -= 1
        u[5, i] += 1
    end
    return nothing
end

# Rates
birth_rate(state, parameters) = parameters.μ * sum(state)

exposure_rate(state, parameters) = parameters.β * state[1] * state[3]
infection_rate(state, parameters) = parameters.a * state[2]
recovery_rate(state, parameters) = parameters.γ * state[3]

exposed_death_rate(state, parameters) = parameters.μ * state[2]
infected_death_rate(state, parameters) = (parameters.α+parameters.μ) * state[3]
recovered_death_rate(state, parameters) = parameters.μ * state[4]


function ∑a(state, parameters, i)
    a_i =  [birth_rate(state, parameters),

            exposure_rate(state, parameters),
            infection_rate(state, parameters),
            recovery_rate(state, parameters),

            exposed_death_rate(state, parameters),
            infected_death_rate(state, parameters),
            recovered_death_rate(state, parameters),
            ]
    # sort!(a_i)
    return sum(a_i[1:i])
end

# μ=0.01, β=0.07, a=0.07, γ=0.006, α = 0.02
μ=0.01; β=0.03; a=0.07; γ=0.006; α=0.02;
const parameters = My_Parameters(μ, β, a, γ, α)

function direct_ssa(initial_condition, parameters)
    u[:, 1] .= initial_condition

    for i in 1:N_SIM
        state = u[:, i]
        u[:, i+1] .= state

        a₀ = ∑a(state, parameters, N_COMPARTMENT)

        r₁ = rand(Float64)
        r₂ = rand(Float64)
        τ  = -log(r₁)/a₀

        # birth
        if r₂*a₀ < ∑a(state, parameters, 1)
            birth!(u, i+1)
        # exposure
        elseif r₂*a₀ < ∑a(state, parameters, 2)
            susceptible_to_exposed!(u, i+1)
        # infection
        elseif r₂*a₀ < ∑a(state, parameters, 3)
            exposed_to_infectious!(u, i+1)
        # recovery
        elseif r₂*a₀ < ∑a(state, parameters, 4)
            infectious_to_resistant!(u, i+1)
        elseif r₂*a₀ < ∑a(state, parameters, 5)
            exposed_death!(u, i+1)
        elseif r₂*a₀ < ∑a(state, parameters, 6)
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
    return plot(t, [S E I R D],
                label=['S' 'E' 'I' 'R' 'D'],
                xlabel = "Time",
                ylabel = "People",
                title="SEIRD Stochastic Simulation",
                linetype=:steppost)
end

const initial_condition = 999, 0, 1, 0, 0
direct_ssa(initial_condition, parameters)
plot_solution(u)
