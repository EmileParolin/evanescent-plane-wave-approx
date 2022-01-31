@testset "PPW" begin
    k = 2
    U = [0, 0, 0, 0, 1]
    result = Dirichlet_sampling(k, U, 20)
    @test result[4] ≈ 0 atol = 1e-12
    @test result[6](0.5,π/4) ≈ 0 atol = 1e-12
end

@testset "EPW" begin
    k = 1
    U = [0, 0, 0, 0, 1]
    @testset "Uniform" begin
        result = Dirichlet_sampling(k, U, 40; smpl_type=uniform_sampling)
        @test result[4] ≈ 0 atol = 1e-8
        @test result[6](0.5,π/4) ≈ 0 atol = 1e-12
    end
    @testset "Sobol" begin
        result = Dirichlet_sampling(k, U, 40; smpl_type=sobol_sampling)
        @test result[4] ≈ 0 atol = 1e-8
        @test result[6](0.5,π/4) ≈ 0 atol = 1e-12
    end
    @testset "Random" begin
        result = Dirichlet_sampling(k, U, 40; smpl_type=random_sampling)
        @test result[4] ≈ 0 atol = 1e-7
        @test result[6](0.5,π/4) ≈ 0 atol = 1e-11
    end
end