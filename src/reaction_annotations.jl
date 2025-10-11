function name_reactions!(model)
    rxn_name_lu = Dict(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "annotations", "reaction_names.csv"),
            types = [String, String],
        ),
    )
    for rid in A.reactions(model)
        if haskey(rxn_name_lu, rid)
            model.reactions[rid].name = rxn_name_lu[rid]
        end
    end
end

function acronym_reactions!(model)
    rxn_acronym_lu =
        Dict(CSV.File(joinpath(pkgdir(@__MODULE__), "data", "annotations", "reaction_acronyms.csv")))
    for rid in A.reactions(model)
        if haskey(rxn_acronym_lu, rid)
            model.reactions[rid].annotations["acronym"] = [rxn_acronym_lu[rid]]
        end
    end
end

function add_reaction_annotations!(model)
    # these are mappings from proteins to annotations (used to justify why a grr is what it is)
    hamap =
        Dict(CSV.File(joinpath(pkgdir(@__MODULE__), "data", "annotations", "hamap.csv")))

    eggnog = CSV.File(joinpath(pkgdir(@__MODULE__), "data", "annotations", "eggnog.csv"))
    eggnog_ec = Dict(
        k => String.(split(v, ",")) for
        (k, v) in zip(eggnog.Gene, eggnog.ECs) if !ismissing(v)
    )


    # these are mappings from rhea ids to annotations (to reaction xrefs)
    xrefs = CSV.File(
        joinpath(pkgdir(@__MODULE__), "data", "annotations", "reaction_xrefs.csv"),
        types = String,
    )
    for row in xrefs
        push!(
            get!(model.reactions[row.Rhea].annotations, String(row.Database), String[]),
            row.Xref,
        )
    end

    qts = Dict( # TODO consider caching
        k => [vs...] for (k, vs) in
        get_quartets([rid for rid in A.reactions(model) if isdigit(first(rid))])
    )

    for rid in A.reactions(model)
        r = model.reactions[rid]
        # assign reaction annos
        if isdigit(first(rid))
            r.annotations["rhea.reaction"] = qts[rid]
        end
        r.annotations["SBO"] = ["SBO_0000176"]

        # assign gene refs
        _grrs = A.reaction_gene_association_dnf(model, rid)
        if !isnothing(_grrs)
            grrs = vcat(_grrs...)
            for gid in grrs
                if haskey(hamap, gid)
                    r.annotations[gid] = [hamap[gid]]
                end
            end
            if any(haskey(eggnog_ec, x) for x in grrs)
                r.annotations["eggnog.ec"] =
                    vcat([eggnog_ec[gid] for gid in grrs if haskey(eggnog_ec, gid)]...)
            end
        end
    end
end
