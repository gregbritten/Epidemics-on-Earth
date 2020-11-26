### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ 18331730-2f51-11eb-1c29-e761eafd08fa
using Pkg

# ╔═╡ d41a64f2-2f54-11eb-0f41-9721260c087c
begin
	using CSV
	using Plots
	using Dates
	using Distributions
	using Printf
	using LaTeXStrings
end

# ╔═╡ 20067550-2f52-11eb-0eca-0d7d597eb358
Pkg.activate("C:/Users/PhytoGreg/.julia/pluto_notebooks")

# ╔═╡ 0bcdcf72-3001-11eb-04d3-d36cbc5366ff
md"""
#### Load the ensemble from the .CSV
"""

# ╔═╡ e6adf2d0-2f54-11eb-275a-a31ffe67c927
begin
	Npop = 10000
	b    = 0.3
	dd = CSV.read(@sprintf "SEIR_ensemble_n=%.0f_b=%.2f.csv" Npop b)
end

# ╔═╡ 241ef040-3001-11eb-225d-c98e9b4be956
md"""
#### Subset the ensemble for a single trial
"""

# ╔═╡ 21263490-2fff-11eb-3f39-dfc8dc9f9575
d = dd[dd.trial .==1,:]

# ╔═╡ e5f562e0-2ff1-11eb-216e-ef112ba183c5
begin
	S = d.X1
	E = d.X2
	I = d.X3
	R = d.X4
end

# ╔═╡ 43a4daae-3001-11eb-1cc8-1b95ff99f3b2
md"""
#### Plot the trial
"""

# ╔═╡ a487a970-2f52-11eb-37d3-85bc423f12b5
begin
	p1 = plot(S)
	p2 = plot(E)
	p3 = plot(I)
	p4 = plot(R)
	plot(p1,p2,p3,p4,labels=["S" "E" "I" "R"])
end

# ╔═╡ 1d3b5db0-2f6e-11eb-0398-417afc82f773
md"""
# SIMULATE TESTING
#### Perfect random sampling
Sampling the population with the probability of a positive test equal to the fraction of infected individuals in the population, according to 
\begin{equation}
	p(t) = \frac{I(t)}{N_{pop}}.
\end{equation}

The simplest case is a constant number of tests performed each day. Below we will investigate time-dependent number of tests performed per day
"""

# ╔═╡ 3705d80e-2f6e-11eb-10de-5fc904afa132
###--CONSTANT NUMBER PER DAY--######
begin
	ptest1 = 0.01 #proportion of the population tested per day
	ptest2 = 0.001 #proportion of the population tested per day
	ntest1 = ptest1*Npop
	ntest2 = ptest2*Npop
	n     = length(S)
	Ntest1 = repeat([convert(Int64,round(ntest1,digits=0))],Int(n)) #time series of 
	Ntest2 = repeat([convert(Int64,round(ntest2,digits=0))],Int(n)) #time series of tests
end

# ╔═╡ cb070eb0-3001-11eb-2697-19d4b838fc93
md"""
#### Plot the time-dependent fraction of the infected individuals
"""

# ╔═╡ dea16800-2f6e-11eb-3a80-51fce93d786f
begin
	Ifrac = I/Npop
	Ifrac[Ifrac.<=0.0] .= 0.0
	plot(Ifrac,size=(500,300),title="Test Probability = I/Npop",legend=false)
end

# ╔═╡ e8be7ab0-3001-11eb-223c-79f147f52638
md"""
#### Plot a simulated time series for positive tests over time
"""

# ╔═╡ a7cb83a0-2f6f-11eb-140a-23ac62fc37a4
begin
	binom1 = Binomial.(Ntest1,Ifrac)
	binom2 = Binomial.(Ntest2,Ifrac)
	t1 = rand.(binom1)
	t2 = rand.(binom2)
	plot(t1,size=(500,300),labels=ntest1,titlefontsize=10, title = @sprintf "Observed positives for different ntest")
	plot!(t2,size=(500,300),labels=ntest2)
end

# ╔═╡ c90d75f0-2f6f-11eb-058b-7ddb7172cc8e
md"""
### Biased sampling
Here we introduce a tendency to sample infected individuals over others. 

This is done by controlling the fraction of non-infected individuals that are available to testing, according to 
\begin{equation}
	p(t) = \frac{I(t)}{I(t) + s(N_{pop}-I(t))}.
\end{equation}
"""

# ╔═╡ bdcf3200-2f6f-11eb-06cf-b90f79cdde88
begin
	s1 = 0.0001
	s2 = 0.01
	psamp1 = I./(I .+ s1.*(Npop.-I))
	psamp2 = I./(I .+ s2.*(Npop.-I))
	#psamp[psamp.<=0.0] .= 0.0
	plot(psamp1,size=(500,300),labels=s1,title="Test Probabilities for Different s")
	plot!(psamp2,labels=s2)
end

# ╔═╡ 632f18a0-3006-11eb-29e7-9b2925e2106e
md"""
#### Simulated time series of observed positive tests for different s
"""

# ╔═╡ b9b0d020-2f6f-11eb-39db-1b51d514d8a9
begin
	binoms1 = Binomial.(Ntest1,psamp1)
	binoms2 = Binomial.(Ntest1,psamp2)
	ts1 = rand.(binoms1)
	ts2 = rand.(binoms2)
	plot(ts1,size=(500,300),labels=s1)
	plot!(ts2,labels=s2)
end

# ╔═╡ a2b7ec50-2f6f-11eb-2074-59bf587112d8
md"""
### TIME-VARYING $N_{test}$
"""

# ╔═╡ d33ed330-2f6e-11eb-3e10-0fd5daee1129
md"""
### Linear increase in testing over time

We can impose time-dependent testing regimes to investigate its role in controlling observed growth rates
"""

# ╔═╡ 97cfc620-3008-11eb-2102-f113a77dcb93
md"""
#### Plot observed positive tests for different Ntest time series
"""

# ╔═╡ d282cd5e-2ff6-11eb-2b3f-55a53ca117fd
begin
	Ntestlin1 = convert.(Int64,round.(LinRange(0.0, 100.0, n),digits=0))
	Ntestlin2 = convert.(Int64,round.(LinRange(0.0, 500.0, n),digits=0))
	Ntestlin3 = convert.(Int64,round.(LinRange(0.0, 1000.0, n),digits=0))
	binomlin1 = Binomial.(Ntestlin1,psamp2)
	binomlin2 = Binomial.(Ntestlin2,psamp2)
	binomlin3 = Binomial.(Ntestlin3,psamp2)
	tlin1 = rand.(binomlin1)
	tlin2 = rand.(binomlin2)
	tlin3 = rand.(binomlin3)
	plin1 = plot(tlin1,size=(500,300))
	plin2 = plot!(tlin2)
	plin3 = plot!(tlin3)
end

# ╔═╡ ebcc64b0-3006-11eb-3d1f-ef9748669139
begin
	plot(Ntestlin1)
	plot!(Ntestlin2)
	plot!(Ntestlin3)
end

# ╔═╡ c954bc12-3007-11eb-1f81-43a8a7504479
md"""
#### Linear increase in s over time
"""

# ╔═╡ d71d2d50-3007-11eb-0832-1964386d1f3a
begin
	s01 = 0.01
	s02 = 0.8
	ss1 = LinRange(0.01, 0.8, n)
	ss2 = LinRange(0.8, 0.8, n)
	psamps1 = I./(I .+ ss1.*(Npop.-I))
	psamps2 = I./(I .+ ss2.*(Npop.-I))
	plot(psamps1,size=(500,300),labels=s01,title="Test Probabilities for Different initial s0 at t0")
	plot!(psamps2,labels=s02)
end

# ╔═╡ c4f065a2-3009-11eb-348a-110c4721dae1
md"""
#### Plot simulated positive test time series for different s0
"""

# ╔═╡ 772e8b30-3009-11eb-02d9-2f0137f50288
begin
	binomss1 = Binomial.(Ntest1,psamps1)
	binomss2 = Binomial.(Ntest1,psamps2)
	tss1 = rand.(binomss1)
	tss2 = rand.(binomss2)
	pss1 = plot(tss1,size=(500,300),label=s01)
	pss2 = plot!(tss2,label=s02)
end

# ╔═╡ 0a379380-300b-11eb-152c-c57d20461642
md"""
## Massachusetts data for example
"""

# ╔═╡ 07fb5f20-300b-11eb-0dfb-b1021b3d1b57
M = CSV.read("d:/google/working/covid/massachusetts-history.csv")

# ╔═╡ 30b2a67e-300b-11eb-39b6-a5aff4ec9b8a
begin
	pos = M.positiveIncrease
	neg = M.negativeIncrease
	tests = pos + neg
	date = M.date
end

# ╔═╡ 3196c4f0-300b-11eb-1e71-3dc51de47008
begin
	const DTM = Union{Date, DateTime} # 2000-01-01 to become 2000.0
	yfrac(dtm::DTM) = (dayofyear(dtm) - 1) / daysinyear(dtm)
	decimaldate(dtm::DTM) = year(dtm) + yfrac(dtm)

	time = zeros(0)
	for i=1:length(date)
		sp = split(date[i],"/")
		dd = Date(parse(Int,sp[3]),parse(Int,sp[1]),parse(Int,sp[2]))
		append!(time,decimaldate(dd))
	end
end

# ╔═╡ 3aaecd80-300b-11eb-1752-3f1f6e9fcf53
begin
	p9 = plot(time,tests)
	tests_N = convert.(Int64, round.((Npop/6.893E6)*tests, digits=0))
	p10 = plot(time,tests_N)
	plot(p9,p10,size=(650,200),labels=["MA tests" "Scaled"])
end

# ╔═╡ Cell order:
# ╠═18331730-2f51-11eb-1c29-e761eafd08fa
# ╠═20067550-2f52-11eb-0eca-0d7d597eb358
# ╠═d41a64f2-2f54-11eb-0f41-9721260c087c
# ╟─0bcdcf72-3001-11eb-04d3-d36cbc5366ff
# ╠═e6adf2d0-2f54-11eb-275a-a31ffe67c927
# ╟─241ef040-3001-11eb-225d-c98e9b4be956
# ╠═21263490-2fff-11eb-3f39-dfc8dc9f9575
# ╠═e5f562e0-2ff1-11eb-216e-ef112ba183c5
# ╟─43a4daae-3001-11eb-1cc8-1b95ff99f3b2
# ╠═a487a970-2f52-11eb-37d3-85bc423f12b5
# ╟─1d3b5db0-2f6e-11eb-0398-417afc82f773
# ╠═3705d80e-2f6e-11eb-10de-5fc904afa132
# ╟─cb070eb0-3001-11eb-2697-19d4b838fc93
# ╠═dea16800-2f6e-11eb-3a80-51fce93d786f
# ╟─e8be7ab0-3001-11eb-223c-79f147f52638
# ╠═a7cb83a0-2f6f-11eb-140a-23ac62fc37a4
# ╟─c90d75f0-2f6f-11eb-058b-7ddb7172cc8e
# ╠═bdcf3200-2f6f-11eb-06cf-b90f79cdde88
# ╟─632f18a0-3006-11eb-29e7-9b2925e2106e
# ╠═b9b0d020-2f6f-11eb-39db-1b51d514d8a9
# ╟─a2b7ec50-2f6f-11eb-2074-59bf587112d8
# ╟─d33ed330-2f6e-11eb-3e10-0fd5daee1129
# ╠═ebcc64b0-3006-11eb-3d1f-ef9748669139
# ╟─97cfc620-3008-11eb-2102-f113a77dcb93
# ╠═d282cd5e-2ff6-11eb-2b3f-55a53ca117fd
# ╟─c954bc12-3007-11eb-1f81-43a8a7504479
# ╠═d71d2d50-3007-11eb-0832-1964386d1f3a
# ╟─c4f065a2-3009-11eb-348a-110c4721dae1
# ╠═772e8b30-3009-11eb-02d9-2f0137f50288
# ╟─0a379380-300b-11eb-152c-c57d20461642
# ╠═07fb5f20-300b-11eb-0dfb-b1021b3d1b57
# ╠═30b2a67e-300b-11eb-39b6-a5aff4ec9b8a
# ╠═3196c4f0-300b-11eb-1e71-3dc51de47008
# ╠═3aaecd80-300b-11eb-1752-3f1f6e9fcf53
