using Distributions
using Plots

############################################
##--PART 3--################################
############################################
k         = [20, 40, 55, 90]
γ         = 1/7
factor    = 100
R_T_MAX   = 12
r_t_range = LinRange(0,R_T_MAX,R_T_MAX*100 + 1)

lam            = zeros(length(r_t_range), length(k) - 1)
likelihood_r_t = zeros(length(r_t_range), length(k) - 1)

for j in 1:(length(k)-1)
    global lam[1:length(r_t_range), j] = k[j] * exp.(γ.*(r_t_range.-1))
end

#poisson = Distributions.Poisson.(k)
poisson = Distributions.Poisson.(lam)

for l in 1:3
    global likelihood_r_t[:,l] = pdf.(poisson[:,l],k[l+1])/sum(pdf.(poisson[:,l],k[l+1]))
end

plot(r_t_range, likelihood_r_t, label=["k = 40" "k = 55" "k = 90"],
     xlabel="R(t)", ylab="Density")
