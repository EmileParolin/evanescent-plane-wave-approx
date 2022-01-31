using Test
using StableApproxEPW

@testset "adaptive-quad-rule.jl" begin
    include("adaptive-quad-rule.jl")
end

@testset "circular-waves.jl" begin
    include("circular-waves.jl")
end

@testset "herglotz-densities.jl" begin
    include("herglotz-densities.jl")
end

@testset "plane-waves.jl" begin
    include("plane-waves.jl")
end

@testset "Dirichlet_sampling.jl" begin
    include("Dirichlet_sampling.jl")
end