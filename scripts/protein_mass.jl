using JSON, AbstractFBCModels
using JSONFBCModels
import AbstractFBCModels as A
using DataFrames, DataFramesMeta, CSV, Statistics
using CairoMakie

model = convert(JSONFBCModels.JSONFBCModel, VibrioNatriegens.build_model())

df = DataFrame(CSV.File("quantitative_proteome.csv"))

gs = A.genes(model)

conds = unique(df.Condition)
mass_tots = Float64[]

for cond in conds

    d = @subset(df, :Condition .== cond)
    d = @combine(groupby(d, :Protein), :mf = mean(:MassFraction))
    d = @subset(d, in.(:Protein, Ref(gs)))
    push!(mass_tots, sum(d.mf))

end

fig = Figure(fontsize=20,)
ax = Axis(fig[1,1], ylabel="Modeled metabolism mass fraction", xticks=(1:6, conds), xlabelfont=:bold,ylabelfont=:bold, xticklabelrotation=-pi/2)
barplot!(ax, 1:6, mass_tots)
hidedecorations!(ax, label=false, ticks=false, ticklabels=false)
fig

