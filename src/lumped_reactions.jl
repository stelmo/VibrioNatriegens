
function biotin_metabolism_lumped_reaction!(model)

    gs = String[]
    ms = String[]

    model.reactions["biotin_lumped1"] = Reaction(;
        name = "Biotin_Lumped_Reaction1",
        stoichiometry = Dict(
            "RHEA-COMP:9955" => -1.0, # Malonyl-[acp] methyl ester
            "RHEA-COMP:78449" => -1.0,  # Malonyl-[acp]
            "CHEBI:15378" => -3.0, # H+
            "CHEBI:16526" => 1.0, # CO2
            "CHEBI:57783" => -1.0, # NADPH
            "CHEBI:58349" => 1.0, # NADP+
            "CHEBI:15377" => 1.0, # H2O
            "CHEBI:57945" => -1.0, # NADH
            "CHEBI:57540" => 1.0, # NAD+
            "dummy1" => 1.0,
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_020332932.1" => 1.0, # R09543, not sure if we can skip it
                    "WP_014232616.1" => 2.0, # R10115
                    "WP_014232490.1" => 2.0, # R10115
                    "WP_020336009.1" => 4.0, # R10116
                    "WP_020334789.1" => 4.0, # R10116
                    "WP_014233979.1" => 4.0, # R10116
                    "WP_014233729.1" => 4.0, # R10116
                    "WP_024372813.1" => 4.0, # R10116
                    "WP_014232754.1" => 4.0, # R10117
                ),
            ),
        ],
        notes = Dict("Source" => ["Lumped reaction"]),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "biotin_lumped1")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "biotin_lumped1")))

    model.reactions["biotin_lumped2"] = Reaction(;
        name = "Biotin_Lumped_Reaction2",
        stoichiometry = Dict(
            "RHEA-COMP:78449" => -1.0,  # Malonyl-[acp]
            "dummy1" => -1.0,
            "CHEBI:15378" => -3.0, # H+
            "CHEBI:16526" => 1.0, # CO2
            "CHEBI:57783" => -1.0, # NADPH
            "CHEBI:58349" => 1.0, # NADP+
            "CHEBI:15377" => 1.0, # H2O
            "CHEBI:57945" => -1.0, # NADH
            "CHEBI:57540" => 1.0, # NAD+
            "RHEA-COMP:10186" => 1.0, # Pimeloyl-[acp]
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_014232616.1" => 2.0, # R10119
                    "WP_014232490.1" => 2.0, # R10119
                    "WP_020336009.1" => 4.0, # R10120
                    "WP_020334789.1" => 4.0, # R10120
                    "WP_014233979.1" => 4.0, # R10120
                    "WP_014233729.1" => 4.0, # R10120
                    "WP_024372813.1" => 4.0, # R10120
                    "WP_014232754.1" => 1.0, # R10121
                ),
            ),
        ],
        notes = Dict("Source" => ["Lumped reaction"]), 
    )
    append!(gs, A.reaction_gene_association_dnf(model, "biotin_lumped2")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "biotin_lumped2")))


    add_genes!(model, gs)
    add_metabolites!(model, ms)
end