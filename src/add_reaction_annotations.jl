
function add_reaction_annotations!(model)
    # these are mappings from proteins to annotations (used to justify why a grr is what it is)
    hamap = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "annotations", "hamap", "hamap.json"))
    eggnog_go = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "annotations", "eggnog", "eggnog_go.json"))
    eggnog_ec = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "annotations", "eggnog", "eggnog_ec.json"))
    ko = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "annotations", "kegg", "ko.json"))
    
    # these are reaction cross references
    ec = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "annotations", "rhea", "ec.json"))
    kegg = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "annotations", "rhea", "kegg.json"))
    bigg = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "annotations", "rhea", "bigg.json"))
    metacyc = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "annotations", "rhea", "metacyc.json"))
    seed = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "annotations", "rhea", "seed.json"))
    reactome = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "annotations", "rhea", "reactome.json"))
    metanetx = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "annotations", "rhea", "metanetx.json"))
    sabiork = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "annotations", "rhea", "sabiork.json"))

    # this needs to get downloaded, since it depends on the current state of the model
    qts = Dict(
        k => [vs...] for (k, vs) in
        get_quartets([rid for rid in A.reactions(model) if isdigit(first(rid))])
    )

    for rid in A.reactions(model)
        r = model.reactions[rid]

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

            if any(haskey(eggnog_go, x) for x in grrs)
                r.annotations["eggnog.go"] =
                    vcat([eggnog_go[gid] for gid in grrs if haskey(eggnog_go, gid)]...)
            end

            if any(haskey(ko, x) for x in grrs)
                r.annotations["kegg.ec"] = vcat([ko[gid] for gid in grrs if haskey(ko, gid)]...)
            end
        end

        if isdigit(first(rid))
            r.annotations["rhea.reaction"] = qts[rid]
        end

        if isdigit(first(rid)) && haskey(metanetx, rid)
            r.annotations["metanetx.reaction"] = metanetx[rid]
        end

        if isdigit(first(rid)) && haskey(seed, rid)
            r.annotations["seed.reaction"] = seed[rid]
        end

        if isdigit(first(rid)) && haskey(sabiork, rid)
            r.annotations["sabiork.reaction"] = sabiork[rid]
        end

        if isdigit(first(rid)) && haskey(reactome, rid)
            r.annotations["reactome.reaction"] = [reactome[rid]]
        end

        if isdigit(first(rid)) && haskey(metacyc, rid)
            r.annotations["metacyc.reaction"] = [metacyc[rid]]
        end

        if haskey(kegg, rid)
            r.annotations["kegg.reaction"] = kegg[rid]
        end

        if haskey(bigg, rid)
            r.annotations["bigg.reaction"] = bigg[rid]
        end

        if haskey(ec, rid)
            r.annotations["rhea.ec"] = ec[rid]
        end

        r.annotations["SBO"] = ["SBO_0000176"]

    end
end
