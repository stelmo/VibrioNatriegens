"""
$(TYPEDSIGNATURES)

Add exchange reactions for basic metabolites
"""
function add_exchanges!(model)
    exrids = JSON.parsefile(joinpath("data", "sources", "ecoli_exchange_rids.json"))

    mids = intersect(collect(values(exrids)), A.metabolites(model))

    for mid in mids
        model.reactions["EX_"*mid] = A.CanonicalModel.Reaction(;
            name = "Exchange $mid",
            stoichiometry = Dict(mid => -1)
        )
    end

end
