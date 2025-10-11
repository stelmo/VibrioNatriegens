using Test, VibrioNatriegens
using AbstractFBCModels
import AbstractFBCModels as A
using HiGHS, COBREXA, ConstraintTrees
import ConstraintTrees as C
using CSV

full_model = VibrioNatriegens.build_model() # load model once

@testset "Vibrio natriegens metabolic model" begin
    include("memote.jl")
    include("issues.jl")
end
