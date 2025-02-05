

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
        if isnothing(rname)
            # option 1
            if !isnothing(grrs)
                ns = String[]
                for grr in grrs
                    rs = [VibrioNatriegens.gene_symbol(model, g) for g in grr if !isnothing(VibrioNatriegens.gene_symbol(model, g))]
                    isempty(rs) && continue
                    u = join(intersect(rs))
                    x = join(filter(!isempty, ([setdiff(r, u) for r in rs])))
                    push!(ns, u * x)
                end
                rname = isempty(ns) ? nothing : join(unique(ns), "-")
            end
            
            # option 2
            if isnothing(rname) && haskey(lu, rid)
                rname = lu[rid]
            end
        end

        model.reactions[rid].name = rname
    end

    # special cases
    model.reactions["54528"].name = "D-ribose 5 phosphate cyclase"
    model.reactions["28659"].name = "Galactosidases"
    model.reactions["22751"].name = "ligK"
    model.genes["WP_269465656.1"].name = "ligK"
    model.reactions["23523"].name = "propionate coa transferase"
    model.reactions["32266"].name = "POP2"
    model.genes["WP_020336055.1"].name = "POP2"
    model.reactions["22491"].name = "gfa"
    model.reactions["20552"].name = "tsdA"
    model.reactions["27488"].name = "allC"
    model.reactions["21371"].name = "pucL"
    model.reactions["10807"].name = "allantoin isomerase"
    model.reactions["26304"].name = "pucL"
    model.reactions["17032"].name = "allB"
    model.reactions["33870"].name = "pucG"
        
    
end
