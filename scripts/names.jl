using VibrioNatriegens
import AbstractFBCModels as A
using CSV, DataFrames, DataFramesMeta


model = VibrioNatriegens.build_model()

df = DataFrame(CSV.File(joinpath("heckman.csv"),))
rids = df.react_id

bigg = JSON.parsefile(joinpath("data", "annotations", "rhea", "bigg.json"))

preferred_name = Dict()
for rid in collect(keys(model.reactions))
    avail = get(bigg, rid, [])
    ns = intersect(avail, rids)
    isempty(ns) && continue
    preferred_name[rid] = first(ns)
end
preferred_name

shortlu = Dict(CSV.File())
for (k, v) in preferred_name
shortlu[k] = v
end

df = DataFrame(Reaction=collect(keys(shortlu)), Name=collect(values(shortlu)))

CSV.write(joinpath( "data", "model", "reaction_shortnames.csv"), df)




