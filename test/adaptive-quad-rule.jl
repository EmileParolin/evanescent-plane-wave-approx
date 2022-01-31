
@testset "ϵ-support" begin
    a, b = ϵsupport(inv, 1)
    @test a ≈ -1 atol = 1e-2
    @test b ≈ +1 atol = 1e-2
end

@testset "Quad Rule" begin
    q = adaptive_G_quad(identity)
    @test q[1] ≈ 0 atol = 1e-14
    q = adaptive_G_quad(identity; a = 0)
    @test q[1] ≈ 1/2 atol = 1e-14
end