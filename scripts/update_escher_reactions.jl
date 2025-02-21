using DataFrames, DataFramesMeta, CSV, JSON, VibrioNatriegens
using AbstractFBCModels
import AbstractFBCModels as A

escher_maps = [
    "amino_acids.json"
    "carbohydrate_metabolism.json"
    "cofactors_etc.json"
    "energy_metabolism.json"
    "lipids.json"
    "nucleotides.json"
]

ems = Dict()
for em in escher_maps
    ems[em] = JSON.parsefile(joinpath("maps", em))
end

em = ems["amino_acids.json"]
em[2]["reactions"]["1"]

model = VibrioNatriegens.build_model()
rxns = Set(A.reactions(model))

for (emid, em) in ems
    emm = em[2]
    for rxn in values(emm["reactions"])
        if rxn["bigg_id"] in rxns
            lb, ub = model.reactions[rxn["bigg_id"]].lower_bound, model.reactions[rxn["bigg_id"]].upper_bound
            if lb < -1 && ub > 1
                rxn["reversibility"] = true
            end
        end
    end

    open(joinpath("maps", "up_"*emid), "w") do io
        JSON.print(io, [em[1], emm])
    end
end
