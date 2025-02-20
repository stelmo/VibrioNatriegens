

function set_default_exchanges!(model)

    default_carbon_source = "CHEBI:15903" # glucose
    
    substrates = [
        "CHEBI:15903" # glucose
        "CHEBI:16189" # so4
        "CHEBI:15379" # o2
        "CHEBI:28938" # nh4(+)
        "CHEBI:43474" # pi
        # "CHEBI:29101" # Na+
    ]

    bidirs = [
        "CHEBI:15377" # H2O
    ]

    for mid in [substrates; bidirs]
        if mid == default_carbon_source
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
    model.reactions["63204"].name = "octadecanoate [acyl-carrier-protein] hydrolase"
    model.reactions["41932"].name = "hexadecanoate [acyl-carrier-protein] hydrolase"
    model.reactions["30123"].name = "tetradecanoate [acyl-carrier-protein] hydrolase"
    model.reactions["30119"].name = "dodecanoate [acyl-carrier-protein] hydrolase"
    model.reactions["30115"].name = "decanoate [acyl-carrier-protein] hydrolase"
    model.reactions["30131"].name = "octanoate [acyl-carrier-protein] hydrolase"
    model.reactions["23592"].name = "yahK"
    model.reactions["19045"].name = "cimA"
    model.reactions["70298"].name = "aspartate aminotransferase"
    model.reactions["18609"].name = "argF"
    model.reactions["33054"].name = "betA"
    model.reactions["45700"].name = "gbcAB"
    model.reactions["45767"].name = "DoeC"    
    model.reactions["19533"].name = "astB"
    model.reactions["16953"].name = "astC"
    model.reactions["19737"].name = "PuuD"
    model.reactions["13636"].name = "PuuA"
    model.reactions["28417"].name = "PuuB"
    model.reactions["42343"].name = "PuuC"
    model.reactions["30710"].name = "spuC"
    model.reactions["15692"].name = "cansdh"
    model.reactions["34118"].name = "nspC"
    model.reactions["13389"].name = "GCDH"
    model.reactions["21384"].name = "PAL"    
    model.reactions["64823"].name = "amidase"
    model.reactions["15449"].name = "HGD"
    model.reactions["16833"].name = "GSR"  
    model.reactions["33791"].name = "sucrose hydrolase"  
    model.reactions["20517"].name = "UDP-N-acetylglucosamine 4-epimerase"
    model.reactions["21036"].name = "acetyl-CoA C-acetyltransferase"
    model.reactions["12921"].name = "5-aminolevulinate synthase"
    model.reactions["54644"].name = "trans-4-hydroxy-L-proline dehydratase"
    model.reactions["21152"].name = "4-hydroxyproline epimerase"
    
end
