
n_people  = Int(1E5)
βs        = range(0.05,0.5,step=0.05)
ntrials   = Int(1E4)

cd("/users/gregorybritten/dropbox/working/covid19/urop/simulations/")

files = readdir()
deleteat!(files, 1)

EX = Array{Float64}(undef, ntrials, length(βs))

for k=1:length(βs)
    β₀ = βs[k]

    str = @sprintf "SEIR_ensemble_n=%.0f_b=%.2f.csv" n_people β₀

    dat = CSV.read(files[1],DataFrame)

    for i=1:10000
        d = filter(:trial => isequal(i),dat)
        I = d[!,"X3"]

        EX[i,k] = findall(x->x==0, I)[1]
    end
end
