using Plots
using CSV
using Distributions

###################PART ONE#####################
#Note --   In .csv file, first column is day count, second row is new cases
#          Make sure that there are titles for both columns
#          If there is a specific file, just make 12-31 into comments and use line 32

country_dict = Dict("AL" => "Alabama",
    "AK" => "Alaska",
    "AR" => "Arkansas",
    "AZ" => "Arizona",
    "BEL" => "Belgium",
    "CA" => "California",
    "CO" => "Colorado",
    "CT" => "Connecticut",
    "DC" => "DC",
    "DE" => "Delaware",
    "FL" => "Florida",
    "FRA" => "France",
    "GA" => "Georgia",
    "HI" => "Hawaii",
    "ID" => "Idaho",
    "IL" => "Illinois",
    "MA" => "Massachusetts",
    "NY" => "New York",
    "SWE" => "Sweden",
    "TX" => "Texas",
    "WA" => "Washington",
    "WI" => "Wisconsin")

println("Abbreviation of the Desired State/Country: \n")
x = readline()

country_name = "$x"

country_data_name = country_dict[country_name]
data_dir = joinpath("..", "data", "nytimes_data", country_data_name * ".csv")

m = CSV.read(data_dir)

#m = CSV.read("./Downloads/COVID Data/$(country_dict[country_name]).csv")
#m = CSV.read("./Downloads/[file_name].csv")

function transpose_data(data)
    days = length(m[:, 1])
    transposed = zeros(Int64, 2, days)
    for row = 1:2
        transposed[row, :] = data[:, row]
    end
    return transposed
end

new_cases = transpose_data(m)

#new_cases = Int64.(zeros(2, length(m[:, 1])))
#for i in 1:2
#    global new_cases[i, :] = m[:, i]
#end

function roll_mean(data::Array, win_size = 7)
    rolling = zeros(length(data), 1)
    for i in 1:length(data)
        if i <= win_size
            rolling[i] = sum(data[1:i])/(i)
            i += 1
        elseif i >= win_size
            rolling[i] = sum(data[(i-win_size):i])/win_size
            i += 1
        end
    end
    return rolling
end

smooth = roll_mean(new_cases[2,:])

plot(1:(length(new_cases[2,:])),
    smooth,
    title = "$country_name New Cases - Smoothed",
    label = "$country_name",
    background_color = :black)

###################PART TWO#####################
#LIKELIHOOD OF R_T WITH A GIVEN K

k = round.(Int, smooth)
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

for l in 1:(length(k)-1)
    global likelihood_r_t[:,l] = pdf.(poisson[:,l],k[l+1]) ./ sum(pdf.(poisson[:,l],k[l+1]))
end

plot(r_t_range,
    likelihood_r_t,
    legend = false,
    ylab = "Probability",
    xlab = "R_t",
    background_color = :black)
title!("Likelihood of R_t given k\n$country_name")

###################PART THREE###################
#BAYESIAN UPDATE TO THE PREVIOUS PART

posteriors = zeros(length(r_t_range), length(k) - 1)
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

plot(r_t_range,
    posteriors,
    legend = false,
    ylab = "Probability",
    xlab = "R_t",
    background_color = :black)
title!("Bayesian Update to Previous Model\n$country_name")

###################PART FOUR####################
#PRODUCES MOST LIKELY R_T VALUE

max_posteriors = zeros((length(k) - 1), 1)
for i in 1:length(max_posteriors)
    global max_posteriors[i] = maximum(posteriors[:, i])
end

ind = zeros((length(k) - 1),1)
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
#33-41
println("\n")


for i in 1:length(ind)
    println("\n The most likely value for R_t for Day $i equals $(r_t_range[Int(ind[i])])")
end


most_likely_values = zeros((length(k) - 1), 1)

for i in 1:length(ind)
    global most_likely_values[i] = r_t_range[Int(ind[i])]
end

plot(1:length(most_likely_values),
    most_likely_values,
    line_color = :orange,
    label = "Most Likely Values")
title!("Graph of R_t over Time\n$country_name")

###################PART FIVE####################
#HIGHEST DENSITY INTERVAL

hdi_index = zeros(Int64, length(ind), 2)
sum_prob = 0
high = .9
low = .1

for l in 1:length(ind)
    for i in 1:length(r_t_range)
        global sum_prob = sum_prob + posteriors[i, l]
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

one_line = zeros(length(most_likely_values), 1)
for i in 1:length(most_likely_values)
    global one_line[i, 1] = 1
end

plot(1:length(most_likely_values),
    hdi_r_t[:, 1],
    line_color = :gray,
    fill = (hdi_r_t[:, 1], hdi_r_t[:, 2], :gray),
    xlab = "Days", ylab = "R_t",
    title = "R_t by Day\n$(country_dict[country_name])",
    label = "HDI", xlim = (5, length(most_likely_values)),
    ylim = (0, 4))
plot!(1:length(most_likely_values),
    most_likely_values,
    line_color = :orange,
    label = "Most Likely Values")
plot!(1:length(most_likely_values),
    hdi_r_t[:, 1],
    line_color = :green,
    label = "Lower")
plot!(1:length(most_likely_values),
    hdi_r_t[:, 2],
    line_color = :purple,
    label = "Upper")
plot!(1:length(most_likely_values),
    one_line,
    label = "R_t = 1")


#Try to make it so k isn't rounded
#Better way of smoothing out new cases line
    #Triangular or Exponential Weighting
#Attempt to add an input element to things
#Add a element that goes to URL and clicks the link

#Story Map for US
    # Start with pictures at each location
    # Next do each link runs a code and displays a plot
    # Then try to expand
#HI and AK don't work well
