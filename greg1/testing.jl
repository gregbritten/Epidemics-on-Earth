
using CSV
using Plots
using Dates
using Distributions
#########################################################
### READ DATA ###########################################
#########################################################
d = CSV.read("D:/google/working/covid/data/SEIR_n2_B2.csv")

####-PICK PARTICULAR ENSEMBLE--#############
#S = d.x5
#E = d.x6
#I = d.x7
#R = d.x8
S = d.x9
E = d.x10
I = d.x11
R = d.x12

p1 = plot(S)
p2 = plot(E)
p3 = plot(I)
p4 = plot(R)
plot(p1,p2,p3,p4)

###--CUMULATIVE TIME SERIES--###############
Ss = 1000.0 .- cumsum(E)
Es = cumsum(E - I)
Is = cumsum(I - R)
Rs = cumsum(R)

p1 = plot(Ss)
p2 = plot(Es)
p3 = plot(Is)
p4 = plot(Rs)

plot(p1,p2,p3,p4)

###################################################################
### TAKE A LOOK AT MASSACHUSETTS DATA #############################
###################################################################
M = CSV.read("d:/google/working/covid/massachusetts-history.csv")

pos = M.positiveIncrease
neg = M.negativeIncrease
tests = pos + neg
date = M.date

###--CONVERT DATE (NEED BETTER WAY TO DO THIS)--###############
const DTM = Union{Date, DateTime} # 2000-01-01 to become 2000.0
yfrac(dtm::DTM) = (dayofyear(dtm) - 1) / daysinyear(dtm)
decimaldate(dtm::DTM) = year(dtm) + yfrac(dtm)

time = zeros(0)
for i=1:length(date)
    sp = split(date[i],"/")
    dd = Date(parse(Int,sp[3]),parse(Int,sp[1]),parse(Int,sp[2]))
    append!(time,decimaldate(dd))
end

p1 = plot(time,tests)
tests_N = convert.(Int64, round.((1000.0/6.893E6)*tests, digits=0))
p2 = plot(time,tests_N)
plot()
#########################################################################
### SIMULATE TESTING ####################################################
#########################################################################

###--CONSTANT NUMBER PER DAY--######
n = length(S)
Ntest = repeat([10],n)

Ifrac = Is/1000.0
plot(Ifrac)
Ifrac[Ifrac.<=0.0] .= 0.0

binom = Binomial.(Ntest,Ifrac)
t1 = rand.(binom)

plot(t1)
