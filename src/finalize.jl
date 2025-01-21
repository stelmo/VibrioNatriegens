

function set_default_exchanges!(model)

    carbon_source = "CHEBI:15903" # glucose
    
    substrates = [
        "CHEBI:16189" # so4
        "CHEBI:15379" # o2
        "CHEBI:28938" # nh4(+)
        "CHEBI:43474" # pi
        "CHEBI:29101" # Na+
    ]

    bidirs = [
        "CHEBI:15377" # H2O
    ]

    for mid in [substrates; bidirs]

        if mid == carbon_source
            lb, ub = (-22.0, 0.0)
        elseif mid in substrates
            lb, ub = (-1000.0, 0.0)
        elseif mid in bidirs
            lb, ub = (-1000.0, 1000.0)
        end

        model.reactions["EX_$mid"].lower_bound = lb
        model.reactions["EX_$mid"].upper_bound = ub
    end

end

function name_reactions!(model)
    df = DataFrame(CSV.File(joinpath("data", "model", "reaction_names.csv")))
    dropmissing!(df)
    lu = Dict(string.(df.RHEA_ID) .=> String.(df.Name))
    for rid in A.reactions(model)
        if isnothing(model.reactions[rid].name) && haskey(lu, rid)
            model.reactions[rid].name = lu[rid]
        end
    end
end
