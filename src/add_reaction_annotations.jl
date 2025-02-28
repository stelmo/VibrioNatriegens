
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
    hamap = Dict(zip(hamap.Protein, hamap.Subunit))

    for rid in A.reactions(model)
        r = model.reactions[rid]
        grrs = A.reaction_gene_association_dnf(model, rid)
        isnothing(grrs) && continue
        grrs = vcat(grrs...)
        for gid in grrs
            r.annotations[gid] = [get(hamap, gid, ""),]
        end
        
    end
end
