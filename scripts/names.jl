using VibrioNatriegens
import AbstractFBCModels as A
using CSV, DataFrames, DataFramesMeta

df = DataFrame(CSV.File(joinpath("data", "annotations", "metanetx", "reac_xref.tsv"), missingstring = ["EMPTY",]))
dropmissing!(df)

df = @combine(
    groupby(df, :meta),
    :xrefs = unique(:ref),
)

bigg = @rtransform(@rsubset(df, startswith(:xrefs, "bigg.reaction")), :xrefs = last(split(:xrefs, ":")))
seed  = @rtransform(@rsubset(df, startswith(:xrefs, "seed.reaction")), :xrefs = last(split(:xrefs, ":")))
metacyc = @rtransform(@rsubset(df, startswith(:xrefs, "metacyc.reaction")), :xrefs = last(split(:xrefs, ":")))
kegg = @rtransform(@rsubset(df, startswith(:xrefs, "kegg.reaction")), :xrefs = last(split(:xrefs, ":")))
sabio = @rtransform(@rsubset(df, startswith(:xrefs, "sabiork.reaction")), :xrefs = last(split(:xrefs, ":")))


rhea = @rtransform(@rsubset(df, startswith(:xrefs, "rhea")), :xrefs = last(split(:xrefs, ":")))
d = Dict(r.xrefs => r.meta for r in eachrow(rhea))

b = Dict()
for r in eachrow(seed)
    push!(get!(b, r.meta, []), r.xrefs)
end

nb = Dict()
for (k, v) in d
    haskey(b, v) || continue
    nb[k] = b[v]
end

open(joinpath("data", "annotations", "rhea", "seed.json"), "w") do io
    JSON.print(io, nb)
end
