using CSV, DataFrames, DataFramesMeta

eggnog = DataFrame(
    CSV.File(joinpath("data", "annotations", "eggnog", "out.emapper.annotations")),
)

@select!(eggnog, :KEGG_ko, :KEGG_Pathway, :KEGG_Module, :KEGG_Reaction, :BRITE, :KEGG_TC)

rids = String[]
for rid in eggnog.KEGG_Reaction
    rid == "-" && continue
    append!(rids, split(rid, ","))
end
unique(rids)

mdf = DataFrame(CSV.File(joinpath("data", "annotations", "grr", "mapped_grrs.csv")))
mrids = unique(mdf.Reaction)

setdiff(mrids, rids)


