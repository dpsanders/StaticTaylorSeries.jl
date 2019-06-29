using StaticTaylorSeries
using Test

@testset "StaticTaylorSeries.jl" begin
    t = StaticTaylor(1, 2, 3)

    @test t + 1 == StaticTaylor(2, 2, 3)
    @test t - 1 == StaticTaylor(0, 2, 3)
    @test t * 2 == StaticTaylor(2, 4, 6)

    @test t^2 == StaticTaylor(1, 4, 10)

    t = StaticTaylor(1.0, 2.0, 3.0)
    @test t / 2 == StaticTaylor(0.5, 1.0, 1.5)
end
