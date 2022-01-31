@testset "Plane Waves" begin
    @test abs(gpw(0,π/8)(0.5, π/4)) ≈ 1 atol = 1e-14
end