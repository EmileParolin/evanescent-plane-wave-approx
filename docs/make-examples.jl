# This file should be executed manually and locally to update the examples

push!(LOAD_PATH,joinpath(@__DIR__, "../src/"))
using Documenter, Literate, StableApproxEPW

# Generating some documentation files with `Literate.jl`
ex_file = joinpath(@__DIR__, "../example/example.jl")
Literate.markdown(ex_file, dirname(ex_file)*"/../docs/src/"; flavor=Literate.DocumenterFlavor())
Literate.notebook(ex_file, dirname(ex_file))