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
    length(k) > 1 && @info("Multiple compound IDs for $mid")

    first(sort(k)) # return smallest index
end

d = Dict(
    exrid => get_kegg_compound_id(exrid) for
    exrid in exchange_rids if !isnothing(get_kegg_compound_id(exrid))
)

open(joinpath("data", "sources", "ecoli_exchange_rids.json"), "w") do io
    JSON.print(io, d)
end



