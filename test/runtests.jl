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

@testset "dirichlet_sampling.jl" begin
    include("dirichlet_sampling.jl")
end