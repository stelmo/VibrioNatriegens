using JSON, AbstractFBCModels
using JSONFBCModels
import AbstractFBCModels as A
using DataFrames, DataFramesMeta, CSV, Statistics

model = A.load("vnat.json")

m = JSON.parsefile(joinpath("maps", escher_maps[1]))
escher_rids = String[]
escher_map_rids = Dict()
for map_id in escher_maps
    m = JSON.parsefile(joinpath("maps", map_id))
    erids = [v["bigg_id"] for v in values(m[2]["reactions"])]
    escher_map_rids[map_id] = erids
    append!(escher_rids, erids)    
end
escher_rids = unique(escher_rids)

df = DataFrame(CSV.File("quantitative_proteome.csv"))
@subset!(df, :Condition .== "Glucose")
df = @combine(groupby(df, :Protein), :mf = mean(:MassFraction))

gs = A.genes(model)
egids = string.(unique(filter(!isnothing, vcat(vcat([A.reaction_gene_association_dnf(model, rid) for rid in escher_rids]...)...))))


# print difference between full model and kegg
kegg = DataFrame(CSV.File("data/annotations/kegg/ko.txt"; header=["Protein", "KO", "Desc"]))
dropmissing!(kegg)
@rsubset!(kegg, !occursin("ribosom", :Desc)) # remove ribosome
ks = kegg.Protein
d = @subset(df, in.(:Protein, Ref(ks)))
sum(d.mf)

mgs = setdiff(ks, gs)

df2 = @subset(df, in.(:Protein, Ref(mgs)))
@orderby(df2, :mf)

df3 = leftjoin(df2, kegg, on=:Protein) # want to see annotations
df3 = @orderby(df3, -:mf)
CSV.write("missing_enzymes.csv", df3)


# print difference between mapped model and kegg
d = @subset(df, in.(:Protein, Ref(egids)))
sum(d.mf)

mgs = setdiff(ks, egids)

df2 = @subset(df, in.(:Protein, Ref(mgs)))
df3 = leftjoin(df2, kegg, on=:Protein)
df3 = @orderby(df3, -:mf)
CSV.write("missing_enzymes_mapped_kegg.csv", df3)

# full model vs mapped model
mgs = setdiff(gs, egids)

df2 = @subset(df, in.(:Protein, Ref(mgs)))
df3 = leftjoin(df2, kegg, on=:Protein)
df3 = @orderby(df3, -:mf)

CSV.write("missing_enzymes_full_mapped.csv", df3)

