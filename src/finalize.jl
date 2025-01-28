

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
        grrs = A.reaction_gene_association_dnf(model, rid)
        rname = A.reaction_name(model, rid)

        # option 1
        if !isnothing(grrs)
            grr = vcat(grrs...)
            rs = [VibrioNatriegens.gene_symbol(model, g) for g in grr if !isnothing(VibrioNatriegens.gene_symbol(model, g))]
            rname = isempty(rs) ? nothing : join(unique(rs))
        end
        
        # option 2
        if isnothing(rname) && haskey(lu, rid)
            rname = lu[rid]
        end

        model.reactions[rid].name = rname
    end
end
