"""
$(TYPEDSIGNATURES)

Build the skeleton of the model automatically using KEGG's BlastKoala output.
"""
function automatic_build()

    kegg_annos = CSV.File(joinpath("data", "annotations", "kegg", "ko.txt"); header=["Protein", "KO", "Annotations"], delim="\t")

    # draft model
    model = A.CanonicalModel.Model()

    # loop through all KEGG's annotations
    for row in kegg_annos
        gid = row.Protein
        # println(gid)
        ismissing(row.KO) && continue

        for k in split(row.KO, "/")
            ko = VibrioNatriegens.get_kegg_orthology(String(k))
            isnothing(ko.reactions) && continue
            
            for r in ko.reactions
                rid = first(split(r, "  "))
                # println(rid)

                kegg_rxn = VibrioNatriegens.get_kegg_reaction(String(rid))
                if haskey(model.reactions, rid)
                    push!(model.reactions[rid].gene_association_dnf, [gid])
                else
                    model.reactions[rid] = A.CanonicalModel.Reaction(;
                        name = kegg_rxn.name,
                        gene_association_dnf = [[gid]],
                        stoichiometry = kegg_rxn.stoichiometry,
                        annotations = Dict(
                            "KEGG_REACTION" => [kegg_rxn.string_stoichiometry],
                        ),
                    )
                end

                if !haskey(model.genes, gid)
                    model.genes[gid] = A.CanonicalModel.Gene(
                        name = gid,
                    )
                end
                               
                for c in keys(kegg_rxn.stoichiometry)
                    cid = c[1:6]
                    kegg_met = VibrioNatriegens.get_kegg_compound(cid)
                    if !haskey(model.metabolites, cid)
                        model.metabolites[cid] = A.CanonicalModel.Metabolite(
                            name = kegg_met.name,
                            formula = JSONFBCModels.parse_formula(kegg_met.formula),
                        )
                    end
                end
            end
        end
    end

    model
end
