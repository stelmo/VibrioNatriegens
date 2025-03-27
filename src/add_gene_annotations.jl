
function add_gene_annotations!(model)
    gene_df = DataFrame(
        CSV.File(joinpath(
            pkgdir(@__MODULE__),
            "data",
            "annotations",
            "genome",
            "gene_annotations.csv",
        )),
    )
    rename_key(k) = begin
        if k == "proteinaccession"
            return "protein.accession"
        elseif k == "preferred_name"
            return "preferred.name"
        else
            return k
        end
    end
    gene_df.SBO = fill("SBO_0000243", length(gene_df.proteinaccession))

    genes_dict = Dict(
        string(first(gdf.proteinaccession)) => Dict(
            rename_key(k) => [String(string(v))] for
            (k, v) in zip(labels(gdf), gdf[1, :]) if !ismissing(v)
        ) for gdf in groupby(gene_df, :proteinaccession)
    )

    # add annotation info
    for gid in A.genes(model)
        g = model.genes[gid]
        g.name = first(genes_dict[gid]["name"])
        g.annotations = genes_dict[gid]
    end

end
