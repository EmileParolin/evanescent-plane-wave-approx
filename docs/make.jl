push!(LOAD_PATH,joinpath(@__DIR__, "../src/"))
using Documenter, Literate, StableApproxEPW

# # Generating some documentation files with `Literate.jl`
# rm("./src/example.md")
# Literate.markdown("../examples/example.jl", "./src/"; flavor=Literate.DocumenterFlavor())

# Generating documentation with `Documenter.jl`
makedocs(
    modules = [StableApproxEPW],
    clean = false,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename="Stable Approx with EPW",
    authors="Emile Parolin",
)

deploydocs(
    repo = "github.com/EmileParolin/evanescent-plane-wave-approx.git",
)
