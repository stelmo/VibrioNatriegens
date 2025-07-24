

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
            lb, ub = (-25.0, 0.0)
        elseif mid in substrates
            lb, ub = (-1000.0, 0.0)
        elseif mid in bidirs
            lb, ub = (-1000.0, 1000.0)
        end

        model.reactions["EX_$mid"].lower_bound = lb
        model.reactions["EX_$mid"].upper_bound = ub
    end

end

function name_genes!(model)
    gene_name_lu = JSON.parsefile(
        joinpath(
            pkgdir(@__MODULE__),
            "data",
            "model",
            "gene_names.json",
        ),
    )
    for gid in A.genes(model)
        if haskey(gene_name_lu, gid)
            model.genes[gid].name = gene_name_lu[gid]
        end
    end
end

function name_reactions!(model)
    df = DataFrame(
        CSV.File(joinpath(pkgdir(@__MODULE__), "data", "model", "reaction_names.csv")),
    )
    dropmissing!(df)
    rid_name = Dict(string.(df.rhea) .=> df.name)
    for rid in A.reactions(model)
        haskey(rid_name, rid) && (model.reactions[rid].name = rid_name[rid])
    end
end

function shortname_reactions!(model)
    shortlu = Dict(CSV.File(joinpath(pkgdir(@__MODULE__), "data", "model", "reaction_shortnames.csv")))
    for (k, v) in shortlu
        model.reactions[k].annotations["shortname"] = [v,]
    end
end
