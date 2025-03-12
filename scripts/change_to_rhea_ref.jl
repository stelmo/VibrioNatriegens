using RheaReactions, CSV, DataFrames, DataFramesMeta

df = DataFrame(CSV.File(joinpath("data", "model", "metabolic_reactions.csv")))

qts = get_quartets(unique(df.RHEA_ID))

lu = Dict(rr => r for (r, v) in qts for rr in v)

qts[]

@rtransform!(df, :RHEA_ID = lu[string(:RHEA_ID)])
CSV.write(joinpath("data", "model", "metabolic_reactions.csv"), df)
