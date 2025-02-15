using Test, VibrioNatriegens
using CSV

using JSON, AbstractFBCModels
using JSONFBCModels
import AbstractFBCModels as A

include("utils.jl")

@testset "Vibrio natriegens metabolic model" begin

    @testset "Mass balances" begin
        for rid in A.reactions(model)
            @test is_mass_balanced(model, rid)
        end
    end

end
