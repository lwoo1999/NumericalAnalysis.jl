using Test

@testset "Integrate" begin
    @testset "romberg" begin
        @test romberg(x -> x^2, -1, 1, 1e-10) ≈ 2/3
        @test romberg(sin, 0, π, 1e-10) ≈ 2
    end

    @testset "int2d" begin
        @test isapprox(int2d((x, y) -> x^2 + y^2, (0, 1, 100), (0, 1, 100)), 2/3; atol=0.01)
        @test isapprox(int2d((x, y) -> sin(x*y), (0, 1, 100), (0, 1, 100)), 0.24; atol=0.01)
    end
end
