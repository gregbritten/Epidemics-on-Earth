using Random
using Plots
################################################
### SIMPLE BIRTH ###############################
################################################
μ    = 0.5
κ    = 0.45
nsim = Int(1E2)
xin  = 1

N0    = 1
N     = N0
tjump = zeros(nsim+1)
state = zeros(nsim+1)
state[1] = N0

for k = 1:nsim
    if N != 0
        w = rand(Float64)
        t  = -log(w)/(μ*N)

        global N += 1
        state[k+1] = N
        tjump[k+1] = t
    else
        k=nsim
    end
end

t = cumsum(tjump)
plot(t,state,linetype=:steppost)
