using AbstractFBCModels, JSONFBCModels, JSON
import AbstractFBCModels as A

model = A.load(joinpath("data", "ecoli", "iML1515.json"))

exchange_rids = last.(split.(filter(startswith("EX_"), A.reactions(model)), "EX_"))

mid = exchange_rids[1]

get_kegg_compound_id(mid) = begin
    k = get(
        model.metabolites[model.metabolite_index[mid]]["annotation"],
        "kegg.compound",
        nothing,
    )
    isnothing(k) && return nothing
    # println(k)
    k
end

d = Dict(
    exrid => get_kegg_compound_id(exrid) for
    exrid in exchange_rids if !isnothing(get_kegg_compound_id(exrid))
)

open(joinpath("e_coli_exchanges.csv"), "w") do io
    write(io, join(["Metabolite", "KeGG"], ";"), "\n")
    for (k, vs) in d
        for v in vs

            v in mids && write(io, join([k, v], ";"), "\n")
        end
    end
end



