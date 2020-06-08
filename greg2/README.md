# Greg2's directory

To use the "environment" in this directory's `Project.toml`, either start julia with

```
$ julia --project
```

or, with julia open, press `]` and then type

```julia
(@v1.4) pkg> activate .
```

Don't forget the "`.`" 

Then type

```julia
julia> include("SIRD_with_differential_equations.jl")
```

to run a script that solves the Susceptible-Infected-Recovered-Deceased model that we've been discussing.
