#use struct to hold parameters
Pkg.add("Plots")
using Plots

N0 = 1e4
S = N0-5
R = 0
D = 0
c = zeros(30)
c[1] = 5
print(c)
print(sum(c))
μ = 0.01
β₀ = 0.2
β = β₀

quar = 90 # quarantine period
socialdist = 0.05 # factor reducting beta during quaranting

tmax = 1000

# c[0:7] exposed c[7:14] infectious, not ill c[14:30] ill, recovering/ deceased
infectious = cat(zeros(6), [0.2, 0.8], ones(30-8), dims=(1,1))
recov = cat(zeros(14), 0.1*ones(14), [0, 0], dims=(1,1))
dec = cat(zeros(14), 0.02*ones(7), 0.1*ones(7), [1,1], dims=(1,1))
# for stochastic model, these would be probabilities, recalculated
# each day with the values above being the expected values [or using better
# numbers]

t = 0
S_r = []
R_r = []
EI_r = []
I_r = []
D_r = []

while t <= tmax
    global c,R,D,β,S,t
    N = sum(c) + S + R
    I = sum(c.*infectious)
    print(N)
    print(I)
    append!(S_r, S)
    append!(R_r, R)
    append!(EI_r, sum(c))
    append!(I_r, I)
    append!(D_r, D)

    newc = β*S/N*I # Number of new exposed [make stochastic]
    R += sum(c.*recov) + c[30]
    D += sum(c.*dec) + c[30]
    @. c += -c*recov - c*dec - μ*c # apply recovery, covid death, natural mortality

    c[2:30] = c[1:29]
    c[1] = newc
    S += -newc+μ*(N0-N)
    t += 1
    # this is the simple social distancing
    if t > 100 && t < 100+quar
        β = socialdist*β₀
    else
        β = β₀
    end
end

t = [0:tmax...]
print(t)
print(length(t)
)
print(EI_r)
print(length(EI_r))
plot(t, EI_r)
