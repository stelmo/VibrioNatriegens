
function add_gene_annotations!(model)

    gene_df = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "annotations",
                "ncbi",
                "refseq_annotations.tsv",
            ),
        ),
    )
    @rename!(gene_df, :ProteinAccession = $"Protein accession", :GeneID = $"Gene ID")
    @select!(gene_df, :Name, :ProteinAccession, :Symbol, :Begin, :End, :Chromosome, :Orientation, :Accession, :GeneID)
    
    # secondary place to lookup gene symbol
    gene_df2 = DataFrame(
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
    @select!(gene_df2, :query, :Preferred_name)
    @rename!(gene_df2, :ProteinAccession = :query)
    leftjoin!(gene_df, gene_df2, on=:ProteinAccession, matchmissing=:notequal)
    def_symbol(s1, s2) = begin
        ismissing(s1) || return s1
        ismissing(s2) && return missing
        s2 == "-" && return missing
        return s2
    end
    @transform!(gene_df, :Symbol = def_symbol.(:Symbol, :Preferred_name))
    @subset!(gene_df, .!ismissing.(:ProteinAccession))

    genes_dict = Dict(
        string(first(gdf.ProteinAccession)) => Dict(
            k => String(string(v)) for (k, v) in zip(labels(gdf), gdf[1,:]) if !ismissing(v)
        ) for gdf in groupby(gene_df, :ProteinAccession)
    )

    molar_masses = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "genome",
                "vnat_molar_masses.csv",
            ),
        )
    )
    mm = Dict(molar_masses.Protein .=> molar_masses.MW)

    # add annotation info
    for gid in A.genes(model)
        g = model.genes[gid]
        g.name = get(genes_dict[gid], "Name", nothing)
        g.symbol = get(genes_dict[gid], "Symbol", nothing)
        g.molarmass = get(mm, gid, nothing)
        g.annotations = Dict(
            "Accession" => [genes_dict[gid]["Accession"],],
            "Chromosome" => [genes_dict[gid]["Chromosome"],],
            "ChromosomeStart" => [genes_dict[gid]["Begin"],],
            "ChromosomeEnd" => [genes_dict[gid]["End"],],
            "Orientation" => [genes_dict[gid]["Orientation"],],
            "ProteinID" => [genes_dict[gid]["ProteinAccession"],],
            "GeneID" => [genes_dict[gid]["GeneID"],],
        )
    end

end
