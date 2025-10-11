
function name_genes!(model)
    gene_name_lu = Dict( # not allowed to have missing here
        CSV.File(joinpath(pkgdir(@__MODULE__), "data", "annotations", "gene_names.csv")),
    )
    for gid in A.genes(model)
        if haskey(gene_name_lu, gid)
            model.genes[gid].name = gene_name_lu[gid]
        end
    end
end

function add_gene_annotations!(model)
    gene_annos = CSV.File(
        joinpath(pkgdir(@__MODULE__), "data", "annotations", "gene_annotations.csv"),
        types = String,
    )
    # add annotation info
    anno_props = propertynames(first(gene_annos))
    for row in gene_annos
        gid = row.LocusTag
        if haskey(model.genes, gid)
            model.genes[gid].annotations = Dict(
                string(k) => [getproperty(row, k)] for
                k in anno_props if !ismissing(getproperty(row, k))
            )
            model.genes[gid].annotations["SBO"] = ["SBO_0000243"]
        end
    end
end
