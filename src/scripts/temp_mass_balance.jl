using VibrioNatriegens, COBREXA, Gurobi
import AbstractFBCModels as A
import ConstraintTrees as C
using DataFrames, CSV, DataFramesMeta, XLSX


model = VibrioNatriegens.build_model(;)

rxns = filter(!startswith("EX_"), A.reactions(model))

unbalanced = String[]
for rid in rxns
    rid = first(rxns)
    rs = A.reaction_stoichiometry(model, rid)
    atoms = Dict()
    for (mid, stoich) in rs
        _mid = first(split(mid, "_"))
        for (a, c) in A.metabolite_formula(model, mid)
            atoms[a] = get(atoms, a, 0) + c * stoich
        end
    end
    atoms
    all(values(atoms) .== 0) || push!(unbalanced, rid)
end
