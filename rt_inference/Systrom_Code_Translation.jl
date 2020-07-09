using Pkg
Pkg.activate("./Downloads/Distributions.jl/")
using Distributions
using Plots

# RUN EACH PART BY ITSELF TO GET THE DESIRED PLOT

###################PART ONE#######################
#POISSON DISTRIBUTION OF CASES

lambdas = [10, 20, 30, 40]
arrays = zeros(71, length(lambdas))
poisson = Distributions.Poisson.(lambdas)
for i in 1:length(lambdas)
    global arrays[1:71, i] = Distributions.pdf(poisson[i], 0:1:70)
end
plot(0:1:70, arrays, xlab = "New Cases", ylab = "Probability", lab = ["10" "20" "30" "40"], background_color = :black)
title!("Poisson Distribution of Cases\n P(k | λ)")

###################PART TWO#######################
#LIKELIHOOD WITH A GIVEN K

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
#LIKELIHOOD OF R_T WITH A GIVEN K

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
    global lam[1:length(r_t_range), j] = k[j] * exp.(γ.*(r_t_range.-1))
end

poisson = Distributions.Poisson.(lam)

for l in 1:3
    global likelihood_r_t[:,l] = pdf.(poisson[:,l],k[l+1])./sum(pdf.(poisson[:,l],k[l+1]))
end

plot(r_t_range, likelihood_r_t, label= ["40" "55" "90"], ylab = "Probability", xlab = "R_t", background_color = :black)
title!("Likelihood of R_t given k")

###################PART FOUR######################
#BAYESIAN UPDATE TO THE PREVIOUS PART

posteriors = zeros(length(r_t_range), length(k) - 1)
for i in 1:(length(k)-1)
    if i == 1
        global posteriors[:, 1] = likelihood_r_t[:,1]
        i += 1
    else
        cycles = 1
        global posteriors[:, i] = likelihood_r_t[:, i]
        while cycles != i
            global posteriors[:, i] = posteriors[:, i] .* likelihood_r_t[:, i - cycles]
            cycles += 1
        end
        global posteriors[:, i] = posteriors[:, i] ./ sum(posteriors[:, i])
        i += 1
    end
end

plot(r_t_range, posteriors, label= ["40" "55" "90"], ylab = "Probability", xlab = "R_t", background_color = :black)
title!("Bayesian Update to Previous Model")

###################PART FIVE######################
#PRODUCES MOST LIKELY R_T VALUE

max_posteriors = zeros(3, 1)
for i in 1:length(max_posteriors)
    global max_posteriors[i] = maximum(posteriors[:, i])
end

ind = zeros(3,1)
for l in 1:(length(max_posteriors))
    for j in 1:(length(r_t_range)-1)
        if posteriors[j, l] != max_posteriors[l]
            j += 1
        elseif posteriors[j, l] == max_posteriors[l]
            global ind[l] = j
            j = 1
            l += 1
            if l <= length(max_posteriors)
                continue
            else
                break
            end
        end
    end
end

println("\n")

for i in 1:length(ind)
    println("\n The most likely value for R_t for Day $i equals $(r_t_range[Int(ind[i])])")
    global p = vline!([r_t_range[Int(ind[i])]])
end
plot(p, label = ["40" "55" "90" "Day 1 Likely" "Day 2 Likely" "Day 3 Likely"])

most_likely_values = zeros(3, 1)

for i in 1:length(ind)
    global most_likely_values[i] = r_t_range[Int(ind[i])]
end

###################PART SIX#######################
#HIGHEST DENSITY INTERVAL

hdi_index = zeros(Int64, length(ind), 2)
sum_prob = 0

for l in 1:length(ind)
    for i in 1:length(r_t_range)
        global sum_prob = sum_prob + posteriors[i, l]
        abs_error = 1e-3
        if .1 - abs_error <= sum_prob <= .1 + abs_error
            global hdi_index[l, 1] = i
            i += 1
        elseif .9 - abs_error <= sum_prob <= .9 + abs_error && l < length(ind)
            global hdi_index[l, 2] = i
            i += 1
        elseif l == 3
            abs_error = 1.6e-3
            if .9 - abs_error <= sum_prob <= .9 + abs_error
                global hdi_index[l, 2] = i
            else
                i += 1
            end
        elseif i == length(r_t_range) && l < length(ind)
            l += 1
            i = 1
            sum_prob = 0
        else
            i += 1
        end
    end
end

hdi_r_t = zeros(length(ind), 2)

for i in 1:length(ind), j in 1:2
    global hdi_r_t[i, j] = r_t_range[hdi_index[i,j]]
end

plot(1:length(most_likely_values), hdi_r_t[:, 1], line_color = :gray, fill = (hdi_r_t[:, 1], hdi_r_t[:, 2], :gray), xlab = "Days", ylab = "R_t", title = "R_t by Day", label = "HDI")
plot!(1:length(most_likely_values), most_likely_values, line_color = :orange, label = "Most Likely Values")
plot!(1:length(most_likely_values), hdi_r_t[:, 1], line_color = :green, label = "Lower")
plot!(1:length(most_likely_values), hdi_r_t[:, 2], line_color = :purple, label = "Upper")
