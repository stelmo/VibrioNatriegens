
function add_electron_transport_chain!(model)

    gs = String[]
    ms = String[]

    model.reactions["nfn"] = Reaction(; # this is speculative
        name = "NAD(P)+ transhydrogenase (ferredoxin)",
        stoichiometry = Dict(
            "CHEBI:33738" => -2.0, # Reduced ferredoxin
            "CHEBI:57945" => -1.0,  # nadh
            "CHEBI:58349" => -2.0, # NADP+
            "CHEBI:15378" => -1.0, # H+
            "CHEBI:33737" => 2.0, # Oxidized ferredoxin
            "CHEBI:57540" => 1.0, # NAD
            "CHEBI:57783" => 2.0, # NADPH
        ),
        lower_bound = -1000.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(; gene_product_stoichiometry = Dict("WP_020336371.1" => 1.0)),
        ], # nfnB, missing nfnA...
        annotations = Dict(
            "KEGG_REACTION" => [
                "2 reduced [2Fe-2S]-[ferredoxin] + NADH + 2 NADP(+) + H(+) = 2 oxidized [2Fe-2S]-[ferredoxin] + NAD(+) + 2 NADPH",
            ],
            "EC" => ["1.6.1.4"],
            "RHEA" => ["RHEA:47000"],
        ),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "nfn")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "nfn")))

    model.reactions["H_ATPsynthase"] = Reaction(;
        name = "F-type H+-transporting ATPase",
        stoichiometry = Dict(
            "CHEBI:30616" => 1, # atp
            "CHEBI:15377" => 1, # water
            "CHEBI:43474" => -1, # pi
            "CHEBI:456216" => -1, # adp
            "CHEBI:15378" => 3.0,  # h+
            "CHEBI:15378_p" => -4.0,  # h+
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_020333565.1" => 1.0,
                    "WP_014233411.1" => 1.0,
                    "WP_014233417.1" => 1.0,
                    "WP_014233414.1" => 3.0,
                    "WP_014233413.1" => 1.0,
                    "WP_014233416.1" => 2.0,
                    "WP_002540812.1" => 10.0,
                    "WP_014233412.1" => 1.0,
                    "WP_014233418.1" => 1.0,
                ),
            ),
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_020333565.1" => 1.0,
                    "WP_020335248.1" => 1.0,
                    "WP_020335246.1" => 1.0,
                    "WP_020335243.1" => 3.0,
                    "WP_020335242.1" => 1.0,
                    "WP_020335244.1" => 2.0,
                    "WP_014234581.1" => 10.0,
                    "WP_020335249.1" => 3.0,
                    "WP_020335247.1" => 1.0,
                ),
            ),
        ],
        annotations = Dict(
            "KEGG_REACTION" =>
                ["ATP + H2O + 4 H(+)(in) = ADP + phosphate + 5 H(+)(out)"],
            "EC" => ["7.1.2.2"],
        ),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "H_ATPsynthase")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "H_ATPsynthase")))

    model.reactions["cyt_bc1"] = Reaction(;
        name = "Cytochrome c oxidase",
        stoichiometry = Dict(
            "CHEBI:24646" => -1.0, # a ubiquinol
            "RHEA_COMP:14399" => -2.0, # Fe(III)-[cytochrome c]
            "CHEBI:15378_p" => 2.0,  # h+ out
            "CHEBI:132124" => 1.0, # a ubiquinone
            "RHEA_COMP:10350" => 2.0, # Fe(II)-[cytochrome c]
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_014230785.1" => 1.0,
                    "WP_014230784.1" => 1.0,
                    "WP_014230783.1" => 1.0,
                ),
            ),
        ],
        annotations = Dict("EC" => ["7.1.1.8"], "RHEA" => ["RHEA:11484"]),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "cyt_bc1")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "cyt_bc1")))

    model.reactions["cyt_c"] = Reaction(;
        name = "Cytochrome c oxidase",
        stoichiometry = Dict(
            "CHEBI:15379" => -1, # o2
            "RHEA_COMP:10350" => -4.0, # Fe(II)-[cytochrome c]
            "CHEBI:15378" => -8.0, # H+
            "RHEA_COMP:14399" => 4.0, # Fe(III)-[cytochrome c]
            "CHEBI:15377" => 2.0, # h2o
            "CHEBI:15378_p" => 4.0, # H+
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_020333285.1" => 1.0,
                    "WP_020333287.1" => 1.0,
                    "WP_020333286.1" => 1.0,
                    "WP_014231840.1" => 1.0,
                ),
            ),
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_020335300.1" => 1.0,
                    "WP_020335301.1" => 1.0,
                    "WP_020335302.1" => 1.0,
                ),
            ),
        ],
        annotations = Dict("EC" => ["7.1.1.9"], "RHEA" => ["RHEA:11436"]),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "cyt_c")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "cyt_c")))

    model.reactions["cyt_c_cbb3"] = Reaction(;
        name = "Cytochrome c oxidase, cbb3-type",
        stoichiometry = Dict(
            "CHEBI:15379" => -1, # o2
            "RHEA_COMP:10350" => -4.0, # Fe(II)-[cytochrome c]
            "CHEBI:15378" => -8.0, # H+
            "RHEA_COMP:14399" => 4.0, # Fe(III)-[cytochrome c]
            "CHEBI:15377" => 2.0, # h2o
            "CHEBI:15378_p" => 4.0, # H+
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_020333285.1" => 1.0,
                    "WP_020333287.1" => 1.0,
                    "WP_020333286.1" => 1.0,
                    "WP_014231840.1" => 1.0,
                ),
            ),
        ],
        annotations = Dict("RHEA" => ["RHEA:11436"], "EC" => ["7.1.1.9"]),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "cyt_c_cbb3")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "cyt_c_cbb3")))

    model.reactions["cyt_bd"] = Reaction(;
        name = "Cytochrome oxidase BD-I",
        stoichiometry = Dict(
            "CHEBI:15377" => 2.0, # h2o
            "CHEBI:15378" => -4.0, # H+
            "CHEBI:15378_p" => 4.0, # H+
            "CHEBI:15379" => -1.0, # o2
            "CHEBI:132124" => 2.0, # a ubiquinone
            "CHEBI:24646" => -2.0, # a ubiquinol
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_000270284.1" => 1.0, # cydX
                    "WP_020332892.1" => 1.0, # cydA
                    "WP_014231356.1" => 1.0, # cydB
                ),
            ),
        ],
        annotations = Dict("RHEA" => ["RHEA:40527"], "EC" => ["7.1.1.7"]),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "cyt_bd")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "cyt_bd")))

    model.reactions["cyt_bo"] = Reaction(;
        name = "Cytochrome oxidase bo3",
        stoichiometry = Dict(
            "CHEBI:24646" => -2.0, # a ubiquinol
            "CHEBI:15379" => -1.0, # o2
            "CHEBI:15378" => -4.0, # H+
            "CHEBI:15377" => 2.0, # h2o
            "CHEBI:15378_p" => 4.0, # H+
            "CHEBI:132124" => 2.0, # a ubiquinone
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_024372973.1" => 1.0, # cyoA
                    "WP_014234376.1" => 1.0, # cyoB
                    "WP_014234375.1" => 1.0, # cyoC
                    "WP_014234374.1" => 1.0, # cyoD
                ),
            ),
        ],
        annotations = Dict("EC" => ["7.1.1.3"], "RHEA" => ["RHEA:30251"]),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "cyt_bo")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "cyt_bo")))

    model.reactions["pnt"] = Reaction(;
        name = "NADPH:NAD+ oxidoreductase H translocase",
        stoichiometry = Dict(
            "CHEBI:57540" => -1.0, # NAD
            "CHEBI:57783" => -1.0, # NADPH
            "CHEBI:15378" => -1.0, # H+
            "CHEBI:57945" => 1.0, # nadh
            "CHEBI:58349" => 1.0, # NADP+
            "CHEBI:15378_p" => 1.0, # H+
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_020334871.1" => 2.0, # pntA
                    "WP_014234219.1" => 2.0, # pntB
                ),
            ),
        ],
        annotations = Dict(
            "KEGG_REACTION" =>
                ["NAD(+) + NADPH + H(+)(in) = NADH + NADP(+) + H(+)(out)"],
            "EC" => ["7.1.1.1"],
            "EXPASY" => ["https://enzyme.expasy.org/EC/7.1.1.1"],
            "RHEA" => ["RHEA:47992"],
        ),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "pnt")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "pnt")))

    add_genes!(model, gs)
    # no need to add metabolites, because they should all already be in the model
    @assert all(in.(ms, Ref(A.metabolites(model))))
end

function add_salt_transducers!(model)
    gs = String[]
    ms = String[]

    model.reactions["Na_ATPsynthase"] = Reaction(;
        name = "F-type Na+-transporting ATPase",
        stoichiometry = Dict(
            "CHEBI:30616" => 1, # atp
            "CHEBI:15377" => 1, # water
            "CHEBI:43474" => -1, # pi
            "CHEBI:456216" => -1, # adp
            "CHEBI:15378" => -1.0,  # h+
            "CHEBI:29101_p" => -4.0, # Na+
            "CHEBI:29101" => 4.0, # Na+
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_020333565.1" => 1.0,
                    "WP_014233411.1" => 1.0,
                    "WP_014233417.1" => 1.0,
                    "WP_014233414.1" => 3.0,
                    "WP_014233413.1" => 1.0,
                    "WP_014233416.1" => 2.0,
                    "WP_002540812.1" => 10.0,
                    "WP_014233412.1" => 1.0,
                    "WP_014233418.1" => 1.0,
                ),
            ),
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_020333565.1" => 1.0,
                    "WP_020335248.1" => 1.0,
                    "WP_020335246.1" => 1.0,
                    "WP_020335243.1" => 3.0,
                    "WP_020335242.1" => 1.0,
                    "WP_020335244.1" => 2.0,
                    "WP_014234581.1" => 10.0,
                    "WP_020335249.1" => 3.0,
                    "WP_020335247.1" => 1.0,
                ),
            ),
        ],
        annotations = Dict(
            "KEGG_REACTION" =>
                ["4 Na(+)(out) + ADP + phosphate + H(+) = 4 Na(+)(in) + ATP + H2O"],
            "EC" => ["7.2.2.1"],
            "RHEA" => ["RHEA:58156"],
        ),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "Na_ATPsynthase")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "Na_ATPsynthase")))

    model.reactions["oad"] = Reaction(;
        name = "oxaloacetate decarboxylase (Na(+) extruding)",
        stoichiometry = Dict(
            "CHEBI:16452" => -1.0, # oxaloacetate
            "CHEBI:29101" => -2.0, # Na+
            "CHEBI:15378" => -1.0, # H+
            "CHEBI:15361" => 1.0, # pyruvate
            "CHEBI:29101_p" => 2.0, # Na+
            "CHEBI:16526" => 1.0, # CO2
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_020333958.1" => 1.0,
                    "WP_020333959.1" => 1.0,
                    "WP_014232970.1" => 1.0,
                ),
            ),
        ],
        annotations = Dict(
            "KEGG_REACTION" =>
                ["oxaloacetate + 2 Na(+)(in) + H(+) = pyruvate + 2 Na(+)(out) + CO2"],
            "EC" => ["7.2.4.2"],
        ),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "oad")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "oad")))


    model.reactions["nqr"] = Reaction(; # like complex 1 but for salt, assume also pumps 4 out
        name = "Na+-transporting NADH:ubiquinone oxidoreductase",
        stoichiometry = Dict(
            "CHEBI:29101_p" => 4.0, # Na+ (n?)
            "CHEBI:57540" => 1.0, # nad
            "CHEBI:24646" => 1.0, # a ubiquinol
            "CHEBI:132124" => -1.0, # a ubiquinone
            "CHEBI:29101" => -4.0, # Na+ (n?)
            "CHEBI:57945" => -1.0, # nadh
            "CHEBI:15378" => -1.0, # h+
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_014232794.1" => 1.0,
                    "WP_014232793.1" => 1.0,
                    "WP_020335166.1" => 1.0,
                    "WP_014232791.1" => 1.0,
                    "WP_014232790.1" => 1.0,
                    "WP_014232789.1" => 1.0,
                ),
            ),
        ],
        annotations = Dict(
            "KEGG_REACTION" => [
                "a ubiquinone + n Na(+)(in) + NADH + H(+) = a ubiquinol + n Na(+)(out) + NAD(+)",
            ],
            "EC" => ["7.2.1.1"],
            "RHEA" => ["RHEA:47748"],
            "EXPASY" => ["https://enzyme.expasy.org/EC/7.2.1.1"],
        ),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "nqr")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "nqr")))

    model.reactions["rnf"] = Reaction(;
        name = "H+/Na+-translocating ferredoxin:NAD+ oxidoreductase",
        stoichiometry = Dict(
            "CHEBI:33738" => -2, # Reduced ferredoxin
            "CHEBI:29101" => -1, # Na+
            "CHEBI:57540" => -1, # NAD
            "CHEBI:15378" => -1, # H+
            "CHEBI:33737" => 2, # Oxidized ferredoxin
            "CHEBI:29101_p" => 1.0, # Na+
            "CHEBI:57945" => 1.0,  # nadh
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_014232537.1" => 1.0,
                    "WP_014232538.1" => 1.0,
                    "WP_020335672.1" => 1.0,
                    "WP_020335673.1" => 1.0,
                    "WP_020335675.1" => 1.0,
                    "WP_020335674.1" => 1.0,
                ),
            ),
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "WP_014234038.1" => 1.0,
                    "WP_020333836.1" => 1.0,
                    "WP_020333837.1" => 1.0,
                    "WP_020333838.1" => 1.0,
                    "WP_020333840.1" => 1.0,
                    "WP_020333839.1" => 1.0,
                ),
            ),
        ],
        annotations = Dict(
            "KEGG_REACTION" => [
                "2 reduced [2Fe-2S]-[ferredoxin] + Na(+)(in) + NAD(+) + H(+) = 2 oxidized [2Fe-2S]-[ferredoxin] + Na(+)(out) + NADH",
            ],
            "EC" => ["7.1.1.11", "7.2.1.2"],
            "EXPASY" => ["https://enzyme.expasy.org/EC/7.2.1.2"],
        ),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "rnf")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "rnf")))

    add_genes!(model, gs)
    # no need to add metabolites, because they should all already be in the model
    @assert all(in.(ms, Ref(A.metabolites(model))))
end

