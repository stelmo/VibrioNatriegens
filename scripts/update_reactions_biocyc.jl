using DataFrames, DataFramesMeta, CSV, RheaReactions
using VibrioNatriegens
import AbstractFBCModels as A



model = VibrioNatriegens.build_model()

biocyc = DataFrame(CSV.File(joinpath("data", "annotations", "rhea", "biocyc_rxns.csv")))
@select!(biocyc, :rheaDir, :metacyc)

rhea_rxn_dir(rxn, qrt) = begin 
    idx = first(indexin([rxn], qrt))
    isnothing(idx) && error("Reaction not found...")
    idx == 1 && return (-1000, 1000)
    idx == 2 && return (0, 1000)
    idx == 3 && return (-1000, 0)
    idx == 4 && return (-1000, 1000)
end

for rid in A.reactions(model)
    qrt = RheaReactions.get_reaction_quartet(parse(Int, rid))
    df = @subset(biocyc, in.(:rheaDir, Ref(rxns)))
    isempty(df) && continue
    lb, ub = rxn_dir(df[1,1], qrt)
    model.reactions[rid].lower_bound = lb
    model.reactions[rid].upper_bound = ub    
end
