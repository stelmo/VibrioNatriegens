
function add_reaction_annotations!(model)
    # pkgdir(@__MODULE__)
    hamap = JSON.parsefile(joinpath("data", "annotations", "hamap", "hamap.json"))
    kegg = JSON.parsefile(joinpath("data", "annotations", "kegg", "kegg.json"))
    eggnog_go = JSON.parsefile(joinpath("data", "annotations", "eggnog", "eggnog_go.json"))
    eggnog_ec = JSON.parsefile(joinpath("data", "annotations", "eggnog", "eggnog_ec.json"))
    ko = JSON.parsefile(joinpath("data", "annotations", "kegg", "ko.json"))
    
    ec = JSON.parsefile(joinpath("data", "annotations", "rhea", "ec.json"))
    

    open(joinpath("data", "annotations", "rhea", "ec.json"), "w") do io
        JSON.print(io, ec)
    end

    ec = DataFrame(CSV.File(joinpath("data", "annotations", "rhea", "ec_rxns.csv")))
    ec = Dict(string(first(gdf.rhea)) => String.(gdf.ec) for gdf in groupby(ec, :rhea))


    qts = Dict(
        k => [vs...] for (k, vs) in
        get_quartets([rid for rid in A.reactions(model) if isdigit(first(rid))])
    )

    metacyc = Dict(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "annotations", "rhea", "biocyc_rxns.csv"),
            drop = [2],
            types = [String, String, String],
        ),
    )

    seed = Dict(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "annotations", "rhea", "seed_rxns.csv"),
            types = [String, String],
        ),
    )

    # kegg_ec_regex = @compile exactly(1, "[EC:")*
    #     between(1,3, DIGIT)*exactly(1, ".")*
    #     between(1,3, DIGIT)*exactly(1, ".")*
    #     between(1,3, DIGIT)*exactly(1, ".")*
    #     between(1,3, DIGIT)*exactly(1, "]")
    kegg_ec_regex =
        r"(?:\[EC:){1}(?:\d){1,3}(?:\.){1}(?:\d){1,3}(?:\.){1}(?:\d){1,3}(?:\.){1}(?:\d){1,3}(?:\]){1}"

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
                _ecs = vcat([ko[gid] for gid in grrs if haskey(ko, gid)]...)
                ecs = [
                    replace(m.match, "[EC:" => "", "]" => "") for ec in _ecs for
                    m in eachmatch(kegg_ec_regex, ec)
                ]
                r.annotations["kegg.ec"] = ecs
            end
        end

        if isdigit(first(rid))
            r.annotations["rhea.reaction"] = qts[rid]
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

        kegg_and_metacyc = [
            get(r.annotations, "kegg.reaction", String[])
            get(r.annotations, "metacyc.reaction", String[])
        ]
        seed_annos = String[]
        for km in kegg_and_metacyc
            haskey(seed, km) && push!(seed_annos, seed[km])
        end
        if !isempty(seed_annos)
            r.annotations["seed.reaction"] = unique(seed_annos)
        end

        r.annotations["SBO"] = ["SBO_0000176"]

    end
end
