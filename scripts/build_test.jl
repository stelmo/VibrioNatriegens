using VibrioNatriegens
using AbstractFBCModels
using JSONFBCModels
import AbstractFBCModels as A
using RheaReactions, DataFrames, DataFramesMeta, CSV

model = VibrioNatriegens.build_model()

VibrioNatriegens.print_reactions(model)
VibrioNatriegens.print_metabolites(model)

m = convert(JSONFBCModels.JSONFBCModel, model)
AbstractFBCModels.save(m,"vnat.json")

# MASS BALANCE
# only the [thioredoxin]-disulfide metabolites give unbalanced reactions
rids = filter(x -> isdigit(first(x)), unique(A.reactions(model)))
unbal_rids = String[]
for rid in rids
    println(rid)
    s = A.reaction_stoichiometry(model, rid)
    m = Dict()
    for (k, v) in s
        for (kk, vv) in A.metabolite_formula(model, k)
            m[kk] = get(m, kk, 0) + vv * v 
        end
    end
    all(values(m) .== 0) || push!(unbal_rids, rid)    
end

df = DataFrame(CSV.File("reactions-model.csv"))
@select!(df, :rid, :Stoichiometry)
@rsubset!(df, occursin("D-glucose",:Stoichiometry))
CSV.write("glucose.csv", df)
