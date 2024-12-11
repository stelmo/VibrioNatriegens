using Test, VibrioNatriegens
using CSV

@testset "Vibrio natriegens metabolic model" begin

    @testset "Utility functions" begin
        include(joinpath("test", "kegg.jl"))
    end

    @testset "Model tests" begin

    end
end
