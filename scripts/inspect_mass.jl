using VibrioNatriegens
using AbstractFBCModels
using JSON, CSV, DataFrames, DataFramesMeta, Statistics
import AbstractFBCModels as A

model = VibrioNatriegens.build_model()

df = DataFrame(CSV.File("quantitative_proteome.csv"))

mf = @combine(
        @groupby(
        @subset(df, :Condition .== "Glucose"),
        :Protein,
    ),
    :mf = mean(:MassFraction) * 1000.0,
)
mf = Dict(zip(mf.Protein, mf.mf))
lu = Dict(gid => get(mf, first(model.genes[gid].annotations["proteinaccession"]), 0.0) for gid in A.genes(model))

d = Dict{String, Float64}()
for rid in A.reactions(model)
    isnothing(A.reaction_gene_association_dnf(model, rid)) && continue
    t = 0.0
    for grr in A.reaction_gene_association_dnf(model, rid)
        for g in grr
            t += get(lu, g, 0)
        end
    end
    d[rid] = t
end

open("vnat_masses.json", "w") do io
    JSON.print(io, d)
end
