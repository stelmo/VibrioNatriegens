

function set_default_exchanges!(model)

    default_carbon_source = "CHEBI:15903" # glucose

    substrates = [
        "CHEBI:15903" # glucose
        "CHEBI:16189" # so4
        "CHEBI:15379" # o2
        "CHEBI:28938" # nh4(+)
        "CHEBI:43474" # pi
    ]

    bidirs = [
        "CHEBI:15377", # H2O
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

function name_reactions_genes!(model)

    for rid in A.reactions(model)
        grrs = A.reaction_gene_association_dnf(model, rid)
        rname = A.reaction_name(model, rid)
        if isnothing(rname)
            # option 1
            if !isnothing(grrs)
                ns = String[]
                for grr in grrs
                    rs = [
                        VibrioNatriegens.gene_symbol(model, g) for
                        g in grr if !isnothing(VibrioNatriegens.gene_symbol(model, g))
                    ]
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
    model.genes["WP_269465656.1"].name = "ligK"
    model.genes["WP_020336055.1"].name = "POP2"

    df = DataFrame(
        CSV.File(joinpath(pkgdir(@__MODULE__), "data", "model", "rename_reactions.csv")),
    )
    rid_name = Dict(df.ID .=> df.Name)
    for rid in A.reactions(model)
        haskey(rid_name, rid) && (model.reactions[rid].name = rid_name[rid])
    end

end

function rename_gene_ids!(model) # from protein id to locus tag
    df = DataFrame(
        CSV.File(joinpath(pkgdir(@__MODULE__), "data", "genome", "locustag_proteinid.csv")),
    )
    lu = Dict(df.ProteinID .=> df.LocusTag)
    ks = collect(keys(model.genes))
    for k in ks
        g = deepcopy(model.genes[k])
        model.genes[lu[k]] = g
        delete!(model.genes, k)
    end

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

metabolite_renamer(s) =
    occursin("CHEBI:", s) ? replace(s, "CHEBI:" => "") : replace(s, ":" => "_") # fallback just replaces :

function rename_metabolite_ids!(model) # from protein id to locus tag
    ks = collect(keys(model.metabolites))

    for k in ks
        occursin(":", k) || continue # not all ids have : in them
        m = deepcopy(model.metabolites[k])
        model.metabolites[metabolite_renamer(k)] = m
        delete!(model.metabolites, k)
    end

    for k in keys(model.reactions)
        st = model.reactions[k].stoichiometry
        kks = collect(keys(st))
        for kk in kks
            occursin(":", kk) || continue # not all ids have : in them
            st[metabolite_renamer(kk)] = st[kk]
            delete!(st, kk)
        end
    end
end

function rename_reaction_ids!(model) # from protein id to locus tag
    ks = collect(keys(model.reactions))
    for k in ks
        if occursin("CHEBI:", k) # not all rids have chebi in them, and they get deleted later on if this control is true
            r = deepcopy(model.reactions[k])
            model.reactions[metabolite_renamer(k)] = r
            delete!(model.reactions, k)
        end
    end
end
