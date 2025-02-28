using VibrioNatriegens
using AbstractFBCModels
import AbstractFBCModels as A
using CSV, DataFrames, DataFramesMeta

model = VibrioNatriegens.build_model()

df = DataFrame(Reaction=String[],Enzyme=String[],Substrates=String[],Products=String[])
for rid in A.reactions(model)
    # println(rid)
    grr = A.reaction_gene_association_dnf(model, rid)
    isnothing(grr) && continue
    if length(grr) == 1
        gid = keys(first(model.reactions[rid].gene_association).gene_product_stoichiometry)
        length(gid) > 1 && continue
        subs = join([first(model.metabolites[k].annotations["InChI"]) for (k,v) in A.reaction_stoichiometry(model, rid) if v < 0],";")
        prods = join([first(model.metabolites[k].annotations["InChI"]) for (k,v) in A.reaction_stoichiometry(model, rid) if v > 0],";")
        push!(df, (rid, model.genes[first(gid)].sequence, subs, subs))
    end
end
df
CSV.write("vnat.csv", df)
# https://deepmolecules.org/Kcat_results/3d0h7g8c1c
