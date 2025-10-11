using GLM, DataFrames, COBREXA, Gurobi, AbstractFBCModels, VibrioNatriegens, ConstraintTrees, DataFramesMeta
import ConstraintTrees as C
import AbstractFBCModels as A
using CairoMakie, AlgebraOfGraphics
using CSV, JSON

model = VibrioNatriegens.build_model()

rxns = []
for (rid, rxn) in model.reactions
    st = abs.(values(rxn.stoichiometry))
    if length(unique(st)) != 1
        push!(rxns, rid)
    end
end

rxns

open("non_canonical_reaction_stoichiometry.txt", "w") do io
    for rid in rxns
        write(io, rid, "\n")
    end
end
