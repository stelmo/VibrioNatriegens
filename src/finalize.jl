
function switch_off_salt_transporters!(model)
    # these reactions cause loops with the H+ transporters
    keep_rxns = [
        "Na_ATPsynthase"
        "oad"
        "nqr"
        "rnf"
        "ANTI_15378_29101_NhaC"
    ]
    for (r, rxn) in model.reactions
        r in keep_rxns && continue
        mets = rxn.stoichiometry
        if "29101" in collect(keys(mets))
            rxn.lower_bound = 0
            rxn.upper_bound = 0
        end
    end
end

function deactivate_antiports!(model)
    #= 
    If both a symport and an antiport reaction exists for the same metabolites,
    then futile cycles are formed.
    Deactivate all antiporters when a corresponding symport is found
    =#
    symps = [join(split(rid,"_")[2:end], "_") for rid in keys(model.reactions) if startswith(rid, "SYM")]
    antis = [join(split(rid,"_")[2:end], "_") for rid in keys(model.reactions) if startswith(rid, "ANTI")]
    for rid in intersect(symps, antis)
        model.reactions["ANTI_$rid"].lower_bound = 0
        model.reactions["ANTI_$rid"].upper_bound = 0
    end
end

function set_default_exchanges!(model)

    default_carbon_source = "15903" # glucose

    substrates = [
        "15903" # glucose
        "16189" # so4
        "15379" # o2
        "43474" # pi
        "29033" # fe(2+)
    ]

    bidirs = [
        "15377" # H2O
        "28938" # nh4(+)
    ]

    for mid in [substrates; bidirs]
        if mid == default_carbon_source
            lb, ub = (-21.0, 0.0)
        elseif mid in substrates
            lb, ub = (-1000.0, 0.0)
        elseif mid in bidirs
            lb, ub = (-1000.0, 1000.0)
        end

        model.reactions["EX_$mid"].lower_bound = lb
        model.reactions["EX_$mid"].upper_bound = ub
    end

end

function check_gene_names(model)
    # must have all gene names assigned
    for (k, v) in model.genes
        isnothing(v.name) && @warn("$k missing gene name")
    end
end

function check_reaction_acronyms(model)
    # must have all acronyms assigned
    for (k, v) in model.reactions
        occursin("_", k) && continue
        isnothing(get(v.annotations, "acronym", nothing)) &&
            @warn("$k missing reaction acronym")
    end
end

function check_rhea_ref(model)
    # model has all 4 rhea reaction ids && id is the reference one
    for (k, v) in model.reactions
        if isdigit(first(k))
            length(v.annotations["rhea.reaction"]) != 4 &&
                @warn("$k does not have all the Rhea references")
            first(sort(v.annotations["rhea.reaction"])) != k &&
                @warn("$k is not the reference rhea reaction")
        end
    end
end
