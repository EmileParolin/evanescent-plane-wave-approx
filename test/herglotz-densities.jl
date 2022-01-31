@testset "Normalization" begin
    k = 2
    @test wz(;k=k)(+1) > 0
    @test wz(;k=k)(-1) > 0
    @test αp(1, logwz(;k=k)) ≈ αp(-1, logwz(;k=k)) atol = 1e-14
    @test αp(4, logwz(;k=k)) ≈ αp(-4, logwz(;k=k)) atol = 1e-14
end

@testset "Herglotz densities" begin
    @test ãp(2)(1,π/2) ≈ 1 / ãp(-2)(1,π/2) atol = 1e-14
end
