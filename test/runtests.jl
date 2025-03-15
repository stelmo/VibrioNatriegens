using Test, VibrioNatriegens
using AbstractFBCModels
import AbstractFBCModels as A
using Gurobi, COBREXA, ConstraintTrees
import ConstraintTrees as C
using DataFrames, CSV, DataFramesMeta

@testset "Vibrio natriegens metabolic model" begin
    include("basics.jl")
end
