push!(LOAD_PATH,joinpath(@__DIR__, "../src/"))
using Documenter, StableApproxEPW

# Generating documentation with `Documenter.jl`
makedocs(
    modules = [StableApproxEPW],
    clean = false,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename="Stable Approx with EPW",
    authors="Emile Parolin",
    pages=[
        "Getting started" => "index.md",
        "Examples" => "example.md"
        "Source" => "implementation.md"
    ],
)

deploydocs(
    repo = "github.com/EmileParolin/evanescent-plane-wave-approx.git",
)
