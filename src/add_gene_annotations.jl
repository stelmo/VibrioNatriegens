
function add_gene_annotations!(model)

    gene_df = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "annotations",
                "ncbi",
                "gene_annotations.csv",
            ),
        ),
    )

    genes_dict = Dict(
        string(first(gdf.ProteinAccession)) => Dict(
            k => String(string(v)) for
            (k, v) in zip(labels(gdf), gdf[1, :]) if !ismissing(v)
        ) for gdf in groupby(gene_df, :ProteinAccession)
    )

    # add annotation info
    for gid in A.genes(model)
        g = model.genes[gid]
        g.name = get(genes_dict[gid], "Name", nothing)
        g.symbol = get(genes_dict[gid], "Symbol", nothing)
        g.molarmass = parse(Float64, genes_dict[gid]["MW"])
        g.sequence = genes_dict[gid]["Sequence"]
        g.annotations = Dict(
            "Accession" => [genes_dict[gid]["Accession"]],
            "Chromosome" => [genes_dict[gid]["Chromosome"]],
            "ChromosomeStart" => [genes_dict[gid]["Begin"]],
            "ChromosomeEnd" => [genes_dict[gid]["End"]],
            "Orientation" => [genes_dict[gid]["Orientation"]],
            "ProteinID" => [genes_dict[gid]["ProteinAccession"]],
            "GeneID" => [genes_dict[gid]["GeneID"]],
            "Symbol" => [get(genes_dict[gid], "Symbol", "")],
        )
    end

end
