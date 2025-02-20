using VibrioNatriegens
using AbstractFBCModels
using JSONFBCModels
using SBMLFBCModels
import AbstractFBCModels as A
using RheaReactions, DataFrames, DataFramesMeta, CSV

model = VibrioNatriegens.build_model()

# MASS BALANCE
rids = filter(x -> isdigit(first(x)), unique(A.reactions(model)))
unbal_rids = String[]
for rid in rids
    s = A.reaction_stoichiometry(model, rid)
    m = Dict()
    for (k, v) in s
        for (kk, vv) in A.metabolite_formula(model, k)
            m[kk] = get(m, kk, 0) + vv * v 
        end
    end
    all(values(m) .== 0) || push!(unbal_rids, rid)    
end
unbal_rids
