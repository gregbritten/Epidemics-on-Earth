using Plots

# Bring over dictionary and constructed functions in other files
include("/Users/morganmayborne/Programming/COMP167/JuliaFiles/Country_dict.jl")
include("/Users/morganmayborne/Programming/COMP167/JuliaFiles/reproduction_number_setup.jl")
include("/Users/morganmayborne/Programming/COMP167/JuliaFiles/reproduction_number_base_functions.jl")

###################INPUTS#####################
#Note --   In .csv file, first column is day count, second row is new cases
#          Make sure that there are titles for both columns
#          If there is a specific file, just make 12-31 into comments and use line 32

# Readlines that require input from the user
input = false
while input == false
    println("\nAbbreviation of the Desired State [2 letter abbreviation] / Country [3 letter country code]:")
    x = uppercase(readline())
    global country_name = "$x"
    try
        attempt = country_dict[country_name]
        global input = true
    catch KeyError
        println("Pick a valid state/country")
        continue
    end
end

input_2 = false
type_list = ["TMA", "E-TMA", "simple", "SMA", "exponential", "EMA", "triangle", "hybrid"]
while input_2 == false
    println("\nType of Moving Average (exponential, simple, E-TMA[recommended], or TMA[recommended]):")
    y = uppercase(readline())
    global mov_avg = "$y"
    if mov_avg in type_list
        global input_2 = true
    else
        println("Please a valid type of moving average")
    end
end

n = 0
try
    global n = country_dict[country_name][3]
catch KeyError
    println("\nn Value (Numerical Value Please)")
    y = readline()
    global n = parse(Int64, y)
end

#####################SETUP#################################

new_cases = transpose_csv(country_name, type = "combined")
posrate = test_posrate(new_cases)

if posrate == true
    cases_data = zeros(2, length(new_cases[2, :]))
    cases_data[1, :] = (new_cases[1, :])[:]
    cases_data[2, :] = n * new_cases[2, :]
    new_cases = cases_data
end
suscept = susceptibles(n, new_cases)
smooth = roll_mean(new_cases, type = "$mov_avg", win_size = 7)
r_t = rollingr_t(n, new_cases, suscept)

plot(1:(length(new_cases[2,:])),
    smooth,
    title = "$country_name New Cases - Smoothed",
    label = "$country_name",
    ylab = "Number of New Cases",
    xlab = "Day Count",
    background_color = :black)

###################PART TWO#####################
#LIKELIHOOD OF R_T WITH A GIVEN K
#k = round.(Int, smooth)
k = smooth
γ = 1/7
factor = 100
R_T_MAX = 12
r_t_range = zeros(R_T_MAX * factor + 1, 1)

for i in 1:length(r_t_range)
    global r_t_range[i, 1] = (R_T_MAX/length(r_t_range))*i
end

likelihood_r_t = base_likelihood(r_t_range, R_T_MAX, k, γ, type = "population")

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
posteriors = bayesian_post(posteriors, likelihood_r_t)

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
max_posteriors = max_post(max_posteriors)

ind = zeros((length(k) - 1),1)
ind = find_max(max_posteriors, ind)

#= For loop that will print the values of each day or days of interest
for i in 1:length(ind)
    println("\n The most likely value for R_t for Day $i equals $(r_t_range[Int(ind[i])])")
end
=#

most_likely_values = zeros((length(k) - 1), 1)
most_likely_values = ind_likely_values(most_likely_values, r_t_range, ind)

plot(1:length(most_likely_values),
    most_likely_values,
    line_color = :orange,
    label = "Most Likely Values")
title!("Graph of R_t over Time\n$country_name")

###################PART FIVE####################
#HIGHEST DENSITY INTERVAL

hdi_index = zeros(Int64, length(ind), 2)

hdi_r_t = high_density_interval(hdi_index)
one_line = horizontal_line(1)


plot(1:length(most_likely_values),
    hdi_r_t[:, 1],
    line_color = :gray,
    fill = (hdi_r_t[:, 1], hdi_r_t[:, 2], :gray),
    xlab = "Days", ylab = "R_t",
    title = "Effective Reproduction Number by Day\n$(country_dict[country_name][1])",
    label = "HDI", xlim = (5, length(most_likely_values)),
    ylim = (0, 4))
plot!(1:length(most_likely_values),
    most_likely_values,
    line_color = :orange,
    label = "Most Likely Values")
#=plot!(1:length(r_t),
    r_t,
    line_color = :black,
    label = "Actual R_t value")=#
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



#Upload into Github
#Check what days the data throws out (RollingFunctions)
#Preferential Testing/Likelihood of Testing Certain People
#Binomial Distribution based on the Likelihood of testing
    #Rate at which Suspectible and Positive People are tested
    #Time Varying
#I/n
#Plot of Percentage of Positive Tests given number of tests known
#Estimate the True number of cases using the percentage of positive Tests
#Scales the raw data based on the proportion of pos. tests
    #Increases (!) OR Decreases (?)

#Story Map for US
    # Start with pictures at each location
    # Next do each link runs a code and displays a plot
    # Then try to expand
#HI, VT, ND, MT, and AK don't work well (Low populations)
