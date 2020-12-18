using Pkg
Pkg.activate("/Users/morganmayborne/Downloads/Distributions.jl/")
using Distributions

function base_likelihood(r_t_range, R_T_MAX, k, γ; type = "population")
    lam = zeros(length(r_t_range), length(k) - 1)
    likelihood_r_t = zeros(length(r_t_range), length(k) - 1)
    if type == "population"
        for j in 1:(length(k)-1)
            global lam[1:length(r_t_range), j] = k[j] * exp.(γ.*(r_t_range.-1))
        end
    elseif type == "proportion"
        for j in 1:(length(k)-1)
            global lam[1:length(r_t_range), j] = k[j]/n * exp.(γ.*(r_t_range.-1))
        end
    end
    lam_sqrt = .√(lam)
    normal = Distributions.Normal.(lam, lam_sqrt)
    poisson = Distributions.Poisson.(lam)

    for l in 1:(length(k)-1)
        global likelihood_r_t[:,l] = pdf.(normal[:,l],k[l+1]) ./ sum(pdf.(normal[:,l],k[l+1]))
    end
    return likelihood_r_t
end

function bayesian_post(posteriors, likelihood_r_t)
    for i in 1:(length(k)-1)
        if i == 1
            global posteriors[:, 1] = likelihood_r_t[:,1]
            i += 1
        elseif i <= 7
            cycles = 1
                global posteriors[:, i] = likelihood_r_t[:, i]
            while cycles != i
                global posteriors[:, i] = posteriors[:, i] .* likelihood_r_t[:, i - cycles]
                cycles += 1
            end
            global posteriors[:, i] = posteriors[:, i] ./ sum(posteriors[:, i])
            i += 1
        else
            cycles = 1
            global posteriors[:, i] = likelihood_r_t[:, i]
            while cycles <= 7
                global posteriors[:, i] = posteriors[:, i] .* likelihood_r_t[:, i - cycles]
                cycles += 1
            end
            global posteriors[:, i] = posteriors[:, i] ./ sum(posteriors[:, i])
            i += 1
        end
    end
    return posteriors
end

function max_post(max_posteriors)
    for i in 1:length(max_posteriors)
        global max_posteriors[i] = maximum(posteriors[:, i])
    end
    return max_posteriors
end


function find_max(max_posteriors, ind)
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
    return ind
end

function ind_likely_values(most_likely_values, r_t_range, ind)
    for i in 1:length(ind)
        global most_likely_values[i] = r_t_range[Int(ind[i])]
    end
    return most_likely_values
end

function high_density_interval(hdi_index; sum_prob = 0, high = .9, low = .1)
    for l in 1:length(ind)
        for i in 1:length(r_t_range)
            sum_prob = sum_prob + posteriors[i, l]
            abs_error = .05
            if low - abs_error <= sum_prob <= low + abs_error
                global hdi_index[l, 1] = i
                i += 1
            elseif high - abs_error <= sum_prob <= high + abs_error
                global hdi_index[l, 2] = i
                i += 1
            elseif i == length(r_t_range)
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
        if hdi_index[i,j] == 0
            global hdi_r_t[i, j] = 0
        else
            global hdi_r_t[i, j] = r_t_range[hdi_index[i,j]]
        end
    end
    return hdi_r_t
end

function horizontal_line(num = 1)
    one_line = zeros(length(most_likely_values), 1)
    for i in 1:length(most_likely_values)
        global one_line[i, 1] = num
    end
    return one_line
end
