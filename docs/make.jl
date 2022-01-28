push!(LOAD_PATH,"../src/")
using Pkg
Pkg.activate("./")
using Documenter, Literate, STrAW

# # Generating some documentation files with `Literate.jl`
# rm("./src/example.md")
# Literate.markdown("../examples/example.jl", "./src/"; flavor=Literate.DocumenterFlavor())

# Generating documentation with `Documenter.jl`
makedocs(sitename="STrAW Documentation.")
