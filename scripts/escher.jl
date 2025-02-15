using JSON, AbstractFBCModels
using JSONFBCModels
import AbstractFBCModels as A

escher_maps = [
    "amino_acids.json"
    "carbohydrate_metabolism.json"
    "cofactors_etc.json"
    "energy_metabolism.json"
    "lipids.json"
    "nucleotides.json"
]

m = JSON.parsefile(joinpath("maps", escher_maps[1]))
m[2]["reactions"]["1"]

escher_rids = String[]
escher_map_rids = Dict()
for map_id in escher_maps
    m = JSON.parsefile(joinpath("maps", map_id))
    erids = [v["bigg_id"] for v in values(m[2]["reactions"])]
    escher_map_rids[map_id] = erids
    append!(escher_rids, erids)    
end
escher_rids = unique(escher_rids)

model = A.load("vnat.json")
rids = filter(x -> isdigit(x[1]), A.reactions(model))

missing_rids = setdiff(rids, escher_rids)

mids = A.metabolites(model)

included_mids = unique(vcat([collect(keys(A.reaction_stoichiometry(model, rid))) for rid in escher_rids]...))

may_as_well_add = String[]
for rid in missing_rids
    mids = collect(keys(A.reaction_stoichiometry(model, rid)))
    (length(mids) - count(in.(mids, Ref(included_mids))) <= 1) && push!(may_as_well_add, rid)
end
may_as_well_add

