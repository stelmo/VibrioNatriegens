using DataFrames, DataFramesMeta, CSV, JSON

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

escher_rids = Int64[]
escher_map_rids = Dict()
for map_id in escher_maps
    m = JSON.parsefile(joinpath("maps", map_id))
    erids = [parse(Int, v["bigg_id"]) for v in values(m[2]["reactions"]) if isdigit(first(v["bigg_id"]))]
    escher_map_rids[map_id] = erids
    append!(escher_rids, erids)    
end
escher_rids = unique(escher_rids)

df = DataFrame(CSV.File(joinpath("data", "model", "all_metabolic_reactions.csv")))
@subset!(df, in.(:RHEA_ID, Ref(escher_rids)))
CSV.write(joinpath("data", "model", "metabolic_reactions.csv"), df)
