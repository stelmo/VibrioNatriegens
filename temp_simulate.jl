using COBREXA, HiGHS, VibrioNatriegens
import AbstractFBCModels as A
import ConstraintTrees as C
using CSV
using DataFrames, DataFramesMeta

model = VibrioNatriegens.build_model()
sol = flux_balance_analysis(model; optimizer = HiGHS.Optimizer)
for (k, v) in sol.fluxes
    startswith(string(k), "EX_") && println(k, " => ", v)
end

rids = collect(keys(model.reactions))

mdf = DataFrame(CSV.File(joinpath("data", "annotations", "grr", "mapped_grrs.csv")))
mrids = unique(mdf.Reaction)

setdiff(rids, mrids)

# write out model
fi = open("model.tab", "w")
for (k, v) in model.reactions
    startswith(k, "EX_") && continue
    ecs =
        haskey(v.annotations, "KEGG_ENZYME") ?
        filter(!isempty, split(first(v.annotations["KEGG_ENZYME"]), " ")) : []
    ms = haskey(v.annotations, "KEGG_MODULE") ? v.annotations["KEGG_MODULE"] : []
    write(fi, join([[k, v.name]; ecs; ms], "\t"), "\n")
    grrs = v.gene_association_dnf
    isnothing(grrs) && continue
    for grr in grrs
        write(fi, "\t", join(grr, "\t"), "\n")
    end
end
close(fi)

using CSV, DataFrames, DataFramesMeta, XLSX, VibrioNatriegens, Gurobi
import COBREXA as X
import AbstractFBCModels as A
import ConstraintTrees as C
