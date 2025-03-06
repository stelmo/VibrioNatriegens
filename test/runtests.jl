using Test, VibrioNatriegens
using AbstractFBCModels
import AbstractFBCModels as A
using Gurobi, COBREXA, ConstraintTrees
import ConstraintTrees as C
using Distributed
using DataFrames, CSV, DataFramesMeta

addprocs(6, exeflags = `--project=$(Base.active_project())`)
@everywhere begin
    using Gurobi, COBREXA, ConstraintTrees, AbstractFBCModels
    import ConstraintTrees as C
    import AbstractFBCModels as A
    using VibrioNatriegens
end

@testset "Vibrio natriegens metabolic model" begin

    include("basics.jl")

end
