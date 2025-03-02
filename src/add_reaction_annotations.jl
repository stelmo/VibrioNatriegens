
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

    bigg = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "rhea",
                "bigg_rxns.csv",
            ),
        ),
    )
    bigg = Dict(
        String.(string(first(gdf.rhea))) =>  String.(gdf.bigg_id) for gdf in groupby(bigg, :rhea)
    )

    kegg = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "rhea",
                "kegg_rxns.csv",
            ),
        ),
    )
    kegg = Dict(
        String.(string(first(gdf.rhea))) =>  String.(gdf.kegg) for gdf in groupby(kegg, :rhea)
    )
    
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
        k => string.(split(v,",")) for (k, v) in zip(String.(eggnog.query), eggnog.EC) if v != "-"
    )
    eggnog_go = Dict(
        k => string.(split(v,",")) for (k, v) in zip(String.(eggnog.query), eggnog.GOs) if v != "-"
    )

    ko = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "annotations",
                "kegg",
                "ko.txt",
            ),
            header=["Protein", "KO", "Desc"],
        ),
    )
    ko = Dict(
        k => v for (k, v) in zip(String.(ko.Protein), ko.Desc) if !ismissing(v)
    )

    ec = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "rhea",
                "ec_rxns.csv",
            ),
        ),
    )
    ec = Dict(
        string(first(gdf.rhea)) => String.(gdf.ec) for gdf in groupby(ec, :rhea)
    )

    protein_locus = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "genome",
                "locustag_proteinid.csv",
            ),
        ),
    )
    protein_locus = Dict(zip(String.(protein_locus.ProteinID),String.(protein_locus.LocusTag)))

    for rid in A.reactions(model)
        r = model.reactions[rid]
        grrs = A.reaction_gene_association_dnf(model, rid)
        isnothing(grrs) && continue
        grrs = vcat(grrs...)

        for gid in grrs
            if haskey(hamap, gid)
                _gid = protein_locus[gid]
                r.annotations[_gid] = [hamap[gid],]
            end
        end
        if haskey(kegg, rid)
            r.annotations["kegg-reaction"] = kegg[rid]
        end
        if haskey(bigg, rid)
            r.annotations["bigg-reaction"] = bigg[rid]
        end
        if haskey(ec, rid)
            r.annotations["rhea-ec"] = ec[rid]
        end
        if any(haskey(eggnog_ec, x) for x in grrs)
            r.annotations["eggnog-ec"] = vcat([eggnog_ec[gid] for gid in grrs if haskey(eggnog_ec, gid)]...)
        end
        if any(haskey(eggnog_go, x) for x in grrs)
            r.annotations["eggnog-go"] = vcat([eggnog_go[gid] for gid in grrs if haskey(eggnog_go, gid)]...)
        end
        if any(haskey(ko, x) for x in grrs)
            r.annotations["kegg-ec"] = vcat([ko[gid] for gid in grrs if haskey(ko, gid)]...)
        end
        
        delete!(r.annotations, "EC")
    end
end
