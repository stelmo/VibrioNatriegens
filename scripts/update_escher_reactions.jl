using DataFrames, DataFramesMeta, CSV, JSON, VibrioNatriegens
using AbstractFBCModels
import AbstractFBCModels as A
using RheaReactions, CSV, DataFrames, DataFramesMeta

model = VibrioNatriegens.build_model()
rids = filter(x -> isdigit(first(x)), A.reactions(model))
qts = get_quartets(rids)
lu = Dict(rr => r for (r, v) in qts for rr in v)

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
em[2]["reactions"]["159"]
ems["nucleotides.json"]

for (emid, em) in ems
    emm = em[2]
    for rxn in values(emm["reactions"])
        rxn["bigg_id"] = get(lu, rxn["bigg_id"], rxn["bigg_id"])
        if occursin("R-", rxn["bigg_id"])
            rxn["bigg_id"] = replace(rxn["bigg_id"], "R-" => "")
        end
        if occursin("-", rxn["bigg_id"])
            rxn["bigg_id"] = replace(rxn["bigg_id"], "-" => "_")
        end

        if rxn["bigg_id"] ∉ A.reactions(model)
            @warn("$(rxn["bigg_id"]) not found in the model!")
            continue
        end
        lb, ub = model.reactions[rxn["bigg_id"]].lower_bound,
        model.reactions[rxn["bigg_id"]].upper_bound
        if lb < -0.1 && ub > 0.1
            rxn["reversibility"] = true
        else
            rxn["reversibility"] = false
        end

        dir = 1.0 # correct the direction 
        if lb < 0 && ub == 0
            dir = -1.0
        end
        st = A.reaction_stoichiometry(model, rxn["bigg_id"])
        rxn_mets = []
        for (k, v) in st
            push!(rxn_mets, Dict("bigg_id" => k, "coefficient" => dir * v))
        end
        rxn["metabolites"] = rxn_mets
    end

    open(joinpath("maps", emid), "w") do io
        JSON.print(io, [em[1], emm])
    end
end
