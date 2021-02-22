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
end

# ╔═╡ 20067550-2f52-11eb-0eca-0d7d597eb358
Pkg.activate("C:/Users/PhytoGreg/.julia/pluto_notebooks")

# ╔═╡ e6adf2d0-2f54-11eb-275a-a31ffe67c927
begin
	d = CSV.read("D:/google/working/covid/data/SEIR_n2_B2.csv")
	Npop = 1E5
end

# ╔═╡ e5f562e0-2ff1-11eb-216e-ef112ba183c5
begin
	S = d.x5
	E = d.x6
	I = d.x7
	R = d.x8
end

# ╔═╡ a487a970-2f52-11eb-37d3-85bc423f12b5
begin
	p1 = plot(S)
	p2 = plot(E)
	p3 = plot(I)
	p4 = plot(R)
	plot(p1,p2,p3,p4)
end

# ╔═╡ af2b8220-2f52-11eb-3ab9-bfd282f43031
begin
	Ss = Npop .- cumsum(E)
	Es = cumsum(E - I)
	Is = cumsum(I - R)
	Rs = cumsum(R)

	p5 = plot(Ss)
	p6 = plot(Es)
	p7 = plot(Is)
	p8 = plot(Rs)

	plot(p5,p6,p7,p8)
end

# ╔═╡ 5f000670-2f6d-11eb-2675-01c72469992d
md"""
# MASSACHUSETTS DATA
"""

# ╔═╡ 6c554380-2f6d-11eb-1c41-09b79df89a08
M = CSV.read("d:/google/working/covid/massachusetts-history.csv")

# ╔═╡ 99c61b50-2f6d-11eb-1872-bf00b1ad1d7a
begin
	pos = M.positiveIncrease
	neg = M.negativeIncrease
	tests = pos + neg
	date = M.date
end

# ╔═╡ b00db07e-2f6d-11eb-3a99-6da8c2cf2d11
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

# ╔═╡ e0f95730-2f6d-11eb-2299-930f35c8d006
begin
	p9 = plot(time,tests)
	tests_N = convert.(Int64, round.((Npop/6.893E6)*tests, digits=0))
	p10 = plot(time,tests_N)
	plot(p9,p10,size=(650,200))
end

# ╔═╡ 1d3b5db0-2f6e-11eb-0398-417afc82f773
md"""
# SIMULATE TESTING
### Perfect random sampling
"""

# ╔═╡ 3705d80e-2f6e-11eb-10de-5fc904afa132
###--CONSTANT NUMBER PER DAY--######
begin
	ptest = 0.01 #percent of the population tested per day
	n     = length(S)
	Ntest = repeat([convert(Int64,round(ptest*Npop,digits=0))],Int(n)) #time series of tests
end

# ╔═╡ 1fb3509e-2f70-11eb-0e08-9d699d23c9c3


# ╔═╡ dea16800-2f6e-11eb-3a80-51fce93d786f
begin
	Ifrac = Is/Npop
	Ifrac[Ifrac.<=0.0] .= 0.0
	plot(Ifrac,size=(500,300))
end

# ╔═╡ c3b42312-2f6f-11eb-1b8c-6dfadb9d5ead


# ╔═╡ a7cb83a0-2f6f-11eb-140a-23ac62fc37a4
begin
	binom = Binomial.(Ntest,Ifrac)
	t1 = rand.(binom)
	plot(t1,size=(500,300))
end

# ╔═╡ c90d75f0-2f6f-11eb-058b-7ddb7172cc8e
md"""
### Biased sampling
"""

# ╔═╡ bdcf3200-2f6f-11eb-06cf-b90f79cdde88
begin
	s = 0.01
	psamp = Is./(Is .+ s.*(Npop.-Is))
	psamp[psamp.<=0.0] .= 0.0
	plot(psamp,size=(500,300))
end

# ╔═╡ b9b0d020-2f6f-11eb-39db-1b51d514d8a9
begin
	binom2 = Binomial.(Ntest,psamp)
	t2 = rand.(binom2)
	plot(t2,size=(500,300))
end

# ╔═╡ a2b7ec50-2f6f-11eb-2074-59bf587112d8
md"""
### TIME-VARYING $N_{test}$
"""

# ╔═╡ 8f200dd0-2f6f-11eb-2999-b9ecec103e38
plot(time,tests_N)

# ╔═╡ 8ae91ae2-2f6f-11eb-31e1-57fb859c35ab
plot(tests_N,pos)

# ╔═╡ d33ed330-2f6e-11eb-3e10-0fd5daee1129


# ╔═╡ Cell order:
# ╠═18331730-2f51-11eb-1c29-e761eafd08fa
# ╠═20067550-2f52-11eb-0eca-0d7d597eb358
# ╠═d41a64f2-2f54-11eb-0f41-9721260c087c
# ╠═e6adf2d0-2f54-11eb-275a-a31ffe67c927
# ╠═e5f562e0-2ff1-11eb-216e-ef112ba183c5
# ╠═a487a970-2f52-11eb-37d3-85bc423f12b5
# ╠═af2b8220-2f52-11eb-3ab9-bfd282f43031
# ╟─5f000670-2f6d-11eb-2675-01c72469992d
# ╠═6c554380-2f6d-11eb-1c41-09b79df89a08
# ╠═99c61b50-2f6d-11eb-1872-bf00b1ad1d7a
# ╠═b00db07e-2f6d-11eb-3a99-6da8c2cf2d11
# ╠═e0f95730-2f6d-11eb-2299-930f35c8d006
# ╟─1d3b5db0-2f6e-11eb-0398-417afc82f773
# ╠═3705d80e-2f6e-11eb-10de-5fc904afa132
# ╠═1fb3509e-2f70-11eb-0e08-9d699d23c9c3
# ╠═dea16800-2f6e-11eb-3a80-51fce93d786f
# ╠═c3b42312-2f6f-11eb-1b8c-6dfadb9d5ead
# ╠═a7cb83a0-2f6f-11eb-140a-23ac62fc37a4
# ╠═c90d75f0-2f6f-11eb-058b-7ddb7172cc8e
# ╠═bdcf3200-2f6f-11eb-06cf-b90f79cdde88
# ╠═b9b0d020-2f6f-11eb-39db-1b51d514d8a9
# ╟─a2b7ec50-2f6f-11eb-2074-59bf587112d8
# ╠═8f200dd0-2f6f-11eb-2999-b9ecec103e38
# ╠═8ae91ae2-2f6f-11eb-31e1-57fb859c35ab
# ╠═d33ed330-2f6e-11eb-3e10-0fd5daee1129
