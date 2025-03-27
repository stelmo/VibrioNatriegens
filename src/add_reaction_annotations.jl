
function add_reaction_annotations!(model)

    hamap = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "annotations",
                "hamap",
                "hamap_subunits.csv",
            ),
        ),
    )
    hamap = Dict(zip(String.(hamap.Protein), hamap.Subunit))

    bigg =
        DataFrame(CSV.File(joinpath(pkgdir(@__MODULE__), "data", "rhea", "bigg_rxns.csv")))
    bigg = Dict(
        String.(string(first(gdf.rhea))) => String.(gdf.bigg_id) for
        gdf in groupby(bigg, :rhea)
    )

    kegg =
        DataFrame(CSV.File(joinpath(pkgdir(@__MODULE__), "data", "rhea", "kegg_rxns.csv")))
    kegg_met = DataFrame(
        CSV.File(joinpath(pkgdir(@__MODULE__), "data", "model", "metabolic_reactions.csv")),
    )
    kegg = Dict(
        String.(string(first(gdf.rhea))) => String.(gdf.kegg) for
        gdf in groupby(kegg, :rhea)
    )
    for (k, v) in zip(string.(kegg_met.RHEA_ID), kegg_met.KEGG_ID)
        if haskey(kegg, k)
            v in kegg[k] || push!(kegg[k], v)
        else
            kegg[k] = [v]
        end
    end

    eggnog = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "annotations",
                "eggnog",
                "out.emapper.annotations",
            ),
        ),
    )
    eggnog_ec = Dict(
        k => string.(split(v, ",")) for
        (k, v) in zip(String.(eggnog.query), eggnog.EC) if v != "-"
    )
    eggnog_go = Dict(
        k => string.(split(v, ",")) for
        (k, v) in zip(String.(eggnog.query), eggnog.GOs) if v != "-"
    )

    reactome = Dict(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "rhea", "reactome_rxns.csv"),
            types = [String, String],
        ),
    )

    ko = DataFrame(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "annotations", "kegg", "ko.txt"),
            header = ["Protein", "KO", "Desc"],
        ),
    )
    ko = Dict(k => v for (k, v) in zip(String.(ko.Protein), ko.Desc) if !ismissing(v))

    ec = DataFrame(CSV.File(joinpath(pkgdir(@__MODULE__), "data", "rhea", "ec_rxns.csv")))
    ec = Dict(string(first(gdf.rhea)) => String.(gdf.ec) for gdf in groupby(ec, :rhea))

    qts = Dict(
        k => [vs...] for (k, vs) in
        get_quartets([rid for rid in A.reactions(model) if isdigit(first(rid))])
    )

    metacyc = Dict(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "rhea", "biocyc_rxns.csv"),
            drop = [2],
            types = [String, String, String],
        ),
    )

    seed = Dict(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__), 
                "data", "rhea", "seed_rxns.csv"),
            types = [String, String],
        ),
    )
    
    kegg_ec_regex = @compile exactly(1, "[EC:")*
        between(1,3, DIGIT)*exactly(1, ".")*
        between(1,3, DIGIT)*exactly(1, ".")*
        between(1,3, DIGIT)*exactly(1, ".")*
        between(1,3, DIGIT)*exactly(1, "]")

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
                ecs = [replace(m.match,"[EC:"=>"","]"=>"") for ec in _ecs for m in eachmatch(kegg_ec_regex,ec)]
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

        kegg_and_metacyc = [get(r.annotations, "kegg.reaction", String[]);get(r.annotations, "metacyc.reaction", String[])]
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
