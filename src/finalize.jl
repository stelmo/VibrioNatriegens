

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
        "28938" # nh4(+)
        "15377" # H2O
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
    # special cases
    gene_name_lu = Dict(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "annotations",
                "genome",
                "gene_names.csv",
            ),
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

function replace_proteinaccession_with_locustag!(model)
    lu = Dict(
        first(model.genes[g].annotations["protein.accession"]) =>
            first(model.genes[g].annotations["locustag"]) for g in A.genes(model)
    )
    # in genes
    for (k, v) in lu
        g = deepcopy(model.genes[k])
        model.genes[v] = g
        delete!(model.genes, k)
    end

    # in reactions
    for k in keys(model.reactions)
        isnothing(model.reactions[k].gene_association) && continue
        for grr in model.reactions[k].gene_association
            kks = collect(keys(grr.gene_product_stoichiometry))
            for kk in kks
                grr.gene_product_stoichiometry[lu[kk]] = grr.gene_product_stoichiometry[kk]
                delete!(grr.gene_product_stoichiometry, kk)
            end
        end
    end
end
