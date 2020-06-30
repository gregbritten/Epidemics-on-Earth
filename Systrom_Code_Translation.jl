using Pkg
Pkg.activate("Distributions.jl/")
using Distributions
using Plots

###################PART ONE#######################

lambdas = [10, 20, 30, 40]
arrays = zeros(71, length(lambdas))
poisson = Distributions.Poisson.(lambdas)
for i in 1:length(lambdas)
    global arrays[1:71, i] = Distributions.pdf(poisson[i], 0:1:70)
end
plot(0:1:70, arrays, xlab = "New Cases", ylab = "Probability", lab = ["10" "20" "30" "40"], background_color = :black)
title!("Poisson Distribution of Cases\n P(k | λ)")

###################PART TWO#######################

k = 20
lam = zeros(90, 1)
for i in 1:length(lam)
    global lam[i, 1] = (45/(length(lam)))*i + .5
end
poisson = Distributions.Poisson(k)
array = zeros(90, 1)
array = Distributions.pdf(poisson, lam)
norm = Distributions.Normal(k, sqrt(k))
array = Distributions.pdf(norm, lam)
plot(1:.5:45.5, array, title= "Likelihood P(k=20 | λ)", xlab = "λ", ylab = "Probability", legend = false, background_color = :black)

###################PART THREE#####################

k = [20, 40, 55, 90]
γ = 1/7
factor = 100
R_T_MAX = 12
r_t_range = zeros(R_T_MAX * factor + 1, 1)
lam = zeros(length(r_t_range), length(k) - 1)
likelihood_r_t = zeros(length(r_t_range), length(k) - 1)

for i in 1:length(r_t_range)
    global r_t_range[i, 1] = (R_T_MAX/length(r_t_range))*i
end

for j in 1:(length(k)-1)
    global lam[1:length(r_t_range), j] = k[j+1] * exp.(γ.*(r_t_range.-1))
end

poisson = Distributions.Poisson.(k)

#=
for l in 1:(length(k)-1)
    global likelihood_r_t[1:length(r_t_range), l] = Distributions.pdf(poisson[l+1], r_t_range)
end
plot(r_t_range, likelihood_r_t, label= "Test")
