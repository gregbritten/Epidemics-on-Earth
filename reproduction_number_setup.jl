using CSV

function transpose_csv(country_name; type = "combined")
    if type == "separate"
        m = CSV.read("/Users/morganmayborne/Programming/COMP167/JuliaFiles/COVID_Data/$(country_dict[country_name][1]).csv")
        new_cases = Int64.(zeros(2, length(m[:, 1])))
        for i in 1:2
            global new_cases[i, :] = m[:, i]
        end
    elseif type == "combined"
        m = CSV.read("/Users/morganmayborne/Programming/COMP167/JuliaFiles/COVID_Data/Data.csv")
        count_m = 0
        #Count missing data
        for i in 1:length(m[:, ((country_dict[country_name])[2]+1)])
            if typeof(m[i, (country_dict[country_name][2]+1)]) == Missing
                count_m += 1
                i += 1
            elseif typeof(m[i, (country_dict[country_name][2]+1)]) != Missing
                i += 1
            end
        end
        # place day count in the new_cases array
        new_cases = zeros(2, length(m[:, (country_dict[country_name][2] + 1)]) - count_m)
        for i in 1:length(new_cases[1, :])
            new_cases[1, i] = i
        end
        # place new case data in the newcases array
        new_cases[2, :] = m[1:length(new_cases[1,:]), (country_dict[country_name][2]+1)]
    end
    return new_cases
end

function test_posrate(new_cases)
    posrate = false
    for i in 1:length(new_cases[2, :])
        if new_cases[2, i] < 1 && new_cases[2, i] != 0
            posrate = true
            break
        end
    end
    return posrate
end

function rollingr_t(n, new_cases, suscept)
    r_0 = 5.7 #recent median r_0 value from July 2020 journal from the CDC
    r_t = zeros(length(new_cases[1, :])-1, 1)
    for j in 1:length(r_t)
        r_t[j] = r_0*suscept[j]/n
    end
    return r_t
end

function susceptibles(n, new_cases)
    suscept = zeros(1, length(new_cases[1, :]))
    for i in 1:length(suscept)
        if i == 1
            suscept[i] = n - new_cases[i]
        else
            suscept[i] = suscept[i-1] - new_cases[i]
        end
    end
    return suscept
end

function roll_mean(data; win_size = 7, type = "simple", α = 0.1)
    data = data[2, :]
    rolling = zeros(length(data))
    if type == "simple" || type == "SMA"
        for i in 1:length(data)
            if i <= win_size
                rolling[i] = sum(data[1:i])/(i)
                i += 1
            elseif i > win_size
                rolling[i] = sum(data[(i-win_size):i])/win_size
                i += 1
            end
        end
    elseif type == "exponential" || type == "EMA"
        weights = zeros(win_size)
        for j in 1:win_size
            global weights[j] = (1 - α)^(j-1)
        end
        weights = reverse(weights)
        for k in 1:length(data)
            if k <= win_size
                rolling[k] = sum(data[1:k] .* weights[1:k])/(k)
                k += 1
            elseif k > win_size
                rolling[k] = sum(data[(k-(win_size-1)):k] .* weights[1:win_size])/win_size
                k += 1
            end
        end
    elseif type == "TMA" || type == "triangle"
        rolling_1 = zeros(length(data))
        for i in 1:length(data)
            if i <= win_size
                rolling_1[i] = sum(data[1:i])/(i)
                i += 1
            elseif i > win_size
                rolling_1[i] = sum(data[(i-win_size):i])/win_size
                i += 1
            end
        end
        for i in 1:length(data)
            if i <= win_size
                rolling[i] = sum(rolling_1[1:i])/(i)
                i += 1
            elseif i > win_size
                rolling[i] = sum(rolling_1[(i-win_size):i])/win_size
                i += 1
            end
        end
    elseif type == "E-TMA" || type == "hybrid"
        weights = zeros(win_size)
        for j in 1:win_size
            global weights[j] = (1 - α)^(j-1)
        end
        weights = reverse(weights)
        rolling_1 = zeros(length(data))
        for i in 1:length(data)
            if i <= win_size
                rolling_1[i] = sum(data[1:i])/(i)
                i += 1
            elseif i > win_size
                rolling_1[i] = sum(data[(i-win_size):i])/win_size
                i += 1
            end
        end
        for k in 1:length(data)
            if k <= win_size
                rolling[k] = sum(rolling_1[1:k] .* weights[1:k])/(k)
                k += 1
            elseif k > win_size
                rolling[k] = sum(rolling_1[(k-(win_size-1)):k] .* weights[1:win_size])/win_size
                k += 1
            end
        end
    end
    return rolling
end

function inputs(input=False)
    if input == False
        println("Abbreviation of the Desired State [2 letter abbreviation] / Country [3 letter country code]:")
        x = uppercase(readline())
        country_name = "$x"

        println("\nType of Moving Average (exponential, simple, E-TMA[recommended], or TMA[recommended]):")
        x = readline()
        mov_avg = "$x"
        return mov_avg, country_name
    end
end
