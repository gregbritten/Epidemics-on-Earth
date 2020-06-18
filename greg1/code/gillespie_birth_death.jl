using Random
using Plots
################################################
### SIMPLE BIRTH-DEATH #########################
################################################
μ    = 0.5
κ    = 0.4
nsim = Int(1E4)
xin  = 1

N0    = 1
N     = N0
tjump = zeros(nsim+1)
state = zeros(nsim+1)
state[1] = N0

for k = 1:nsim
    if N != 0
        w1 = rand(Float64)
        t  = -log(w1)/((μ+κ)*N)
        w2 = rand(Float64)
        if w2 < μ/(μ+κ)
            global N += 1
        else
            global N -= 1
        end
        state[k+1] = N
        tjump[k+1] = t
    else
        k=nsim
    end
end

t = cumsum(tjump)
plot(t,state,linetype=:steppost)
