using Distances,
      Pathogen,
      Random,
      Plots,
      Plots.PlotMeasures;

# Generate population
n = 100
risks = DataFrame(x = rand(Uniform(0, 15), n),
                  y = rand(Uniform(0, 30), n),
                  riskfactor1 = rand(Gamma(), n))

# Precalculate distances
dists = [euclidean([risks[i, :x];
                    risks[i, :y]],
                   [risks[j, :x];
                    risks[j, :y]]) for i = 1:n, j = 1:n]

pop = Population(risks, dists)

# Define risk functions/TN-ILM structure
function _constant(θ::Vector{Float64}, pop::Population, i::Int64)
  return θ[1]
end

function _one(θ::Vector{Float64}, pop::Population, i::Int64)
  return 1.0
end

function _linear(θ::Vector{Float64}, pop::Population, i::Int64)
  return θ[1] * pop.risks[i, :riskfactor1]
end

function _powerlaw(θ::Vector{Float64}, pop::Population, i::Int64, k::Int64)
  d = pop.distances[k, i]
  return d^(-θ[1])
end

rf = RiskFunctions{SEIR}(_constant, # sparks function
                        _one, # susceptibility function
                        _powerlaw, # infectivity function
                        _one, # transmissability function
                        _linear,
                        _linear) # removal function

# Parametrize risk functions for simulation
rparams = RiskParameters{SEIR}([0.0001], # sparks function parameter(s)
                              Float64[], # susceptibility function parameter(s)
                              [4.0], # infectivity function parameter(s)
                              Float64[], # transmissibility function parameter(s)
                              [0.1],
                              [0.1]) # removal function parameter(s)

starting_states = [State_I; fill(State_S, n-1)]

# Initialize Simulation
sim = Simulation(pop, starting_states, rf, rparams)

# Simulate!
simulate!(sim, tmax=200.0)

# p1 = plot(sim.events, 0.0, 200.0, legendfont=font(6), xaxis=font(10), bottom_margin=30px)
p1 = plot(sim.events)
