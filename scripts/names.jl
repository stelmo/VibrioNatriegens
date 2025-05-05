using VibrioNatriegens
import AbstractFBCModels as A
using CSV, DataFrames, DataFramesMeta

df = DataFrame(CSV.File(joinpath("data", "metanetx", "reac_xref.tsv"), missingstring = ["EMPTY",]))
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

model = VibrioNatriegens.build_model()
model.reactions["13237"].annotations

rhea_metanet = Dict()
for (rid, rxn) in model.reactions
    a = rxn.annotations
    # rhea_metanet[rid] = 
end
