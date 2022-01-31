@testset "Normalization" begin
    k = 2
    @test βp(1; k=k) ≈ βp(-1; k=k) atol = 1e-14
    @test βp(4; k=k) ≈ βp(-4; k=k) atol = 1e-14
end

@testset "Circular waves" begin
    k = 2
    @test b̃p(2; k=k)(1,π/2) ≈ +conj(b̃p(-2; k=k)(1,π/2)) atol = 1e-14
    @test b̃p(3; k=k)(1,π/2) ≈ -conj(b̃p(-3; k=k)(1,π/2)) atol = 1e-14
    @test bp(2; k=k)(1,π/2) ≈ +conj(bp(-2; k=k)(1,π/2)) atol = 1e-14
    @test bp(3; k=k)(1,π/2) ≈ -conj(bp(-3; k=k)(1,π/2)) atol = 1e-14
end