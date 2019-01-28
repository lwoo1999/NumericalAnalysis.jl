using Test
using LinearAlgebra

A = [3  2  5  4  6;
     2  1  3 -7  8;
     5  3  2  5 -4;
     4 -7  5  1  3;
     6  8 -4  3  8]

@testset "Eigen" begin
    @testset "eigenValues" begin
        evs = eigenValues(A, 1e-10)
        @test isapprox(sort(evs), [-12.7606, -4.06195, 4.57495, 11.0234, 16.2242]; atol = 1e-4)
    end

    @testset "eigenVectors" begin
        evs = eigenVectors(A, 1e-10)
        normalize!.(evs)
        for i in eachindex(evs)
            if evs[i][1] < 0
                evs[i] = -evs[i]
            end
        end
        sort!(evs, by=first)
        
        exact = [
            [0.0892401, 0.591435, -0.44077, 0.525459, -0.414554], 
            [0.103969, 0.546299, 0.607112, -0.422398, -0.379135], 
            [0.411867, -0.342988, 0.542907, 0.615544, -0.197695], 
            [0.469878, 0.430087, 0.0775824, 0.101, 0.760276], 
            [0.768639, -0.221691, -0.369284, -0.395446, -0.25942]
        ]
        @test all(xy -> isapprox(xy...; atol = 1e-4), zip(evs, exact))
    end
end