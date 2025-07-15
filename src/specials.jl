
function add_electron_transport_chain!(model)

    gs = String[]
    ms = String[]

    model.reactions["nfn"] = Reaction(; # this is speculative - block 
        name = "NAD(P)+ transhydrogenase (ferredoxin)",
        stoichiometry = Dict(
            "33738" => -2.0, # Reduced ferredoxin
            "57945" => -1.0,  # nadh
            "58349" => -2.0, # NADP+
            "15378" => -1.0, # H+
            "33737" => 2.0, # Oxidized ferredoxin
            "57540" => 1.0, # NAD
            "57783" => 2.0, # NADPH
        ),
        lower_bound = 0.0,
        upper_bound = 0.0,
        gene_association = [
            X.Isozyme(; gene_product_stoichiometry = Dict("PN96_RS05215" => 1.0)),
        ], # nfnB, missing nfnA...
        annotations = Dict(
            "rhea-reaction-description" => [
                "2 reduced [2Fe-2S]-[ferredoxin] + NADH + 2 NADP(+) + H(+) = 2 oxidized [2Fe-2S]-[ferredoxin] + NAD(+) + 2 NADPH",
            ],
            "rhea.ec" => ["1.6.1.4"],
            "rhea.reaction" => ["47000"],
        ),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "nfn")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "nfn")))

    model.reactions["H_ATPsynthase"] = Reaction(;
        name = "F-type H+-transporting ATPase",
        stoichiometry = Dict(
            "30616" => 1, # atp
            "15377" => 1, # water
            "43474" => -1, # pi
            "456216" => -1, # adp
            "15378" => 3.0,  # h+
            "15378_p" => -4.0,  # h+
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "PN96_RS13585" => 1.0,
                    "PN96_RS13605" => 1.0,
                    "PN96_RS13570" => 1.0,
                    "PN96_RS13590" => 3.0,
                    "PN96_RS13595" => 1.0,
                    "PN96_RS13580" => 2.0,
                    "PN96_RS13575" => 10.0,
                    "PN96_RS13600" => 1.0,
                    "PN96_RS13565" => 1.0,
                ),
            ),
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "PN96_RS13585" => 1.0,
                    "PN96_RS22635" => 1.0,
                    "PN96_RS22625" => 1.0,
                    "PN96_RS22610" => 3.0,
                    "PN96_RS22605" => 1.0,
                    "PN96_RS22615" => 2.0,
                    "PN96_RS22620" => 10.0,
                    "PN96_RS22640" => 3.0,
                    "PN96_RS22630" => 1.0,
                ),
            ),
        ],
        annotations = Dict(
            "rhea-reaction-description" =>
                ["ATP + H2O + 4 H(+)(in) = ADP + phosphate + 5 H(+)(out)"],
            "rhea.ec" => ["7.1.2.2"],
            "rhea.reaction" => ["57720"],
        ),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "H_ATPsynthase")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "H_ATPsynthase")))

    model.reactions["cyt_bc1"] = Reaction(;
        name = "Cytochrome c oxidase",
        stoichiometry = Dict(
            "24646" => -1.0, # a ubiquinol
            "29034" => -2.0, # Fe(III)-[cytochrome c]
            "15378_p" => 2.0,  # h+ out
            "132124" => 1.0, # a ubiquinone
            "29033" => 2.0, # Fe(II)-[cytochrome c]
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "PN96_RS11440" => 1.0, # petC
                    "PN96_RS11445" => 1.0, # petB
                    "PN96_RS11450" => 1.0, # petA
                ),
            ),
        ],
        annotations = Dict("rhea.ec" => ["7.1.1.8"], "rhea.reaction" => ["11484"]),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "cyt_bc1")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "cyt_bc1")))

    model.reactions["cyt_c"] = Reaction(;
        name = "Cytochrome c oxidase",
        stoichiometry = Dict(
            "15379" => -1, # o2
            "29033" => -4.0, # Fe(II)-[cytochrome c]
            "15378" => -8.0, # H+
            "29034" => 4.0, # Fe(III)-[cytochrome c]
            "15377" => 2.0, # h2o
            "15378_p" => 4.0, # H+
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "PN96_RS22940" => 1.0, # coxC
                    "PN96_RS22950" => 1.0, # coxA
                    "PN96_RS22955" => 1.0, # coxB
                    "PN96_RS22945" => 1.0, # coxG
                    "PN96_RS22920" => 1.0, # cox15
                ),
            ),
            X.Isozyme(; # Cytochrome c oxidase, cbb3-type
                gene_product_stoichiometry = Dict(
                    "PN96_RS05815" => 1.0, # ccoN
                    "PN96_RS05820" => 1.0, # ccoO
                    "PN96_RS05825" => 1.0, # ccoQ
                    "PN96_RS05830" => 1.0, # ccoP
                    "PN96_RS05840" => 1.0, # ccoI
                    "PN96_RS05845" => 1.0, # ccoS
                ),
            ),
        ],
        annotations = Dict("rhea.ec" => ["7.1.1.9"], "rhea.reaction" => ["11436"]),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "cyt_c")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "cyt_c")))

    model.reactions["cyt_bd"] = Reaction(; # does not pump protons https://doi.org/10.1016/j.bbabio.2011.06.016
        name = "Cytochrome oxidase BD-I",
        stoichiometry = Dict(
            "15377" => 2.0, # h2o
            "15378" => -4.0, # H+
            "15378_p" => 4.0, # H+
            "15379" => -1.0, # o2
            "132124" => 2.0, # a ubiquinone
            "24646" => -2.0, # a ubiquinol
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "PN96_RS08265" => 1.0, # cydX
                    "PN96_RS08270" => 1.0, # cydA
                    "PN96_RS08275" => 1.0, # cydB
                    "PN96_RS07430" => 1.0, # cydD
                    "PN96_RS07435" => 1.0, # cydC
                ),
            ),
        ],
        annotations = Dict("rhea.reaction" => ["40527"], "rhea.ec" => ["7.1.1.7"]),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "cyt_bd")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "cyt_bd")))

    model.reactions["cyt_bo"] = Reaction(; # https://www.pnas.org/doi/10.1073/pnas.2106750118
        name = "Cytochrome oxidase bo3",
        stoichiometry = Dict(
            "24646" => -2.0, # a ubiquinol
            "15379" => -1.0, # o2
            "15378" => -4.0, # H+
            "15377" => 2.0, # h2o
            "15378_p" => 4.0, # H+
            "132124" => 2.0, # a ubiquinone
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "PN96_RS21725" => 1.0, # cyoA
                    "PN96_RS21720" => 1.0, # cyoB
                    "PN96_RS21715" => 1.0, # cyoC
                    "PN96_RS21710" => 1.0, # cyoD
                    "PN96_RS22915" => 1.0, # cyoE
                ),
            ),
        ],
        annotations = Dict("rhea.ec" => ["7.1.1.3"], "rhea.reaction" => ["30251"]),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "cyt_bo")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "cyt_bo")))

    model.reactions["pnt"] = Reaction(;
        name = "NADPH:NAD+ oxidoreductase H translocase",
        stoichiometry = Dict(
            "57540" => -1.0, # NAD
            "57783" => -1.0, # NADPH
            "15378" => -1.0, # H+
            "57945" => 1.0, # nadh
            "58349" => 1.0, # NADP+
            "15378_p" => 1.0, # H+
        ),
        lower_bound = -1000.0,
        upper_bound = 0.0, # https://doi.org/10.1016/j.jmb.2005.07.022
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "PN96_RS20820" => 2.0, # pntA
                    "PN96_RS20825" => 2.0, # pntB
                ),
            ),
        ],
        annotations = Dict(
            "rhea-reaction-description" =>
                ["NAD(+) + NADPH + H(+)(in) = NADH + NADP(+) + H(+)(out)"],
            "rhea.ec" => ["7.1.1.1"],
            "rhea.reaction" => ["47992"],
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
            "30616" => 1, # atp
            "15377" => 1, # water
            "43474" => -1, # pi
            "456216" => -1, # adp
            "15378" => -1.0,  # h+
            "29101_p" => -4.0, # Na+
            "29101" => 4.0, # Na+
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "PN96_RS13585" => 1.0,
                    "PN96_RS13605" => 1.0,
                    "PN96_RS13570" => 1.0,
                    "PN96_RS13590" => 3.0,
                    "PN96_RS13595" => 1.0,
                    "PN96_RS13580" => 2.0,
                    "PN96_RS13575" => 10.0,
                    "PN96_RS13600" => 1.0,
                    "PN96_RS13565" => 1.0,
                ),
            ),
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "PN96_RS13585" => 1.0,
                    "PN96_RS22635" => 1.0,
                    "PN96_RS22625" => 1.0,
                    "PN96_RS22610" => 3.0,
                    "PN96_RS22605" => 1.0,
                    "PN96_RS22615" => 2.0,
                    "PN96_RS22620" => 10.0,
                    "PN96_RS22640" => 3.0,
                    "PN96_RS22630" => 1.0,
                ),
            ),
        ],
        annotations = Dict(
            "rhea-reaction-description" =>
                ["4 Na(+)(out) + ADP + phosphate + H(+) = 4 Na(+)(in) + ATP + H2O"],
            "rhea.ec" => ["7.2.2.1"],
            "rhea.reaction" => ["58156"],
        ),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "Na_ATPsynthase")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "Na_ATPsynthase")))

    model.reactions["oad"] = Reaction(; 
        name = "oxaloacetate decarboxylase (Na(+) extruding)",
        stoichiometry = Dict(
            "16452" => -1.0, # oxaloacetate
            "29101" => -2.0, # Na+
            "15378" => -1.0, # H+
            "15361" => 1.0, # pyruvate
            "29101_p" => 2.0, # Na+
            "16526" => 1.0, # CO2
        ),
        lower_bound = -1000.0, # the delta D for oac <-> pyr + co2 is -32 kJ/mol, hence this is likely pumping Na
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "PN96_RS01245" => 1.0, # oadB
                    "PN96_RS01240" => 1.0, # oadA
                    "PN96_RS01235" => 1.0, # oadG
                ),
            ),
        ],
        annotations = Dict(
            "rhea-reaction-description" =>
                ["oxaloacetate + 2 Na(+)(in) + H(+) = pyruvate + 2 Na(+)(out) + CO2"],
            "rhea.ec" => ["7.2.4.2"],
            "rhea.reaction" => ["57724"],
        ),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "oad")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "oad")))


    model.reactions["nqr"] = Reaction(; # like complex 1 but for salt, assume also pumps 4 out
        name = "Na+-transporting NADH:ubiquinone oxidoreductase",
        stoichiometry = Dict(
            "29101_p" => 4.0, # Na+ (n?)
            "57540" => 1.0, # nad
            "24646" => 1.0, # a ubiquinol
            "132124" => -1.0, # a ubiquinone
            "29101" => -4.0, # Na+ (n?)
            "57945" => -1.0, # nadh
            "15378" => -1.0, # h+
        ),
        lower_bound = 0.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "PN96_RS02145" => 1.0,
                    "PN96_RS02150" => 1.0,
                    "PN96_RS02155" => 1.0,
                    "PN96_RS02160" => 1.0,
                    "PN96_RS02165" => 1.0,
                    "PN96_RS02170" => 1.0,
                ),
            ),
        ],
        annotations = Dict(
            "rhea-reaction-description" => [
                "a ubiquinone + n Na(+)(in) + NADH + H(+) = a ubiquinol + n Na(+)(out) + NAD(+)",
            ],
            "rhea.ec" => ["7.2.1.1"],
            "rhea.reaction" => ["47748"],
        ),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "nqr")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "nqr")))

    model.reactions["rnf"] = Reaction(;
        name = "H+/Na+-translocating ferredoxin:NAD+ oxidoreductase",
        stoichiometry = Dict(
            "33738" => -2, # Reduced ferredoxin
            "29101" => -1, # Na+
            "57540" => -1, # NAD
            "15378" => -1, # H+
            "33737" => 2, # Oxidized ferredoxin
            "29101_p" => 1.0, # Na+
            "57945" => 1.0,  # nadh
        ),
        lower_bound = -1000.0,
        upper_bound = 1000.0,
        gene_association = [
            X.Isozyme(;
                gene_product_stoichiometry = Dict(
                    "PN96_RS03320" => 1.0, # rnfA
                    "PN96_RS03315" => 1.0, # rnfB
                    "PN96_RS03310" => 1.0, # rnfC
                    "PN96_RS03305" => 1.0, # rnfD
                    "PN96_RS03295" => 1.0, # rnfE
                    "PN96_RS03300" => 1.0, # rnfG
                ),
            ),
            X.Isozyme(; # these are not well annotated
                gene_product_stoichiometry = Dict(
                    "PN96_RS19620" => 1.0,
                    "PN96_RS19625" => 1.0,
                    "PN96_RS19630" => 1.0,
                    "PN96_RS19635" => 1.0,
                    "PN96_RS19640" => 1.0,
                    "PN96_RS19645" => 1.0,
                    "PN96_RS19650" => 1.0,
                ),
            ),
        ],
        annotations = Dict(
            "rhea-reaction-description" => [
                "2 reduced [2Fe-2S]-[ferredoxin] + Na(+)(in) + NAD(+) + H(+) = 2 oxidized [2Fe-2S]-[ferredoxin] + Na(+)(out) + NADH",
            ],
            "rhea.ec" => ["7.2.1.2"],
            "rhea.reaction" => ["46800"],
        ),
    )
    append!(gs, A.reaction_gene_association_dnf(model, "rnf")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "rnf")))

    # block most of the transporters, unsure about stoichiometry and they cause loops
    model.reactions["ANTI_15378_29101_NhaA"] = Reaction(
        name = "Antiporter Na/H (NhaA)",
        stoichiometry = Dict(
            "15378" => -2.0, # H+
            "29101_p" => -1.0, # Na+
            "29101" => 1.0, # Na+
            "15378_p" => 2.0, # H+
        ),
        objective_coefficient = 0.0,
        lower_bound = 0,
        upper_bound = 0,
        gene_association = [
            X.Isozyme(; gene_product_stoichiometry = Dict("PN96_RS07530" .=> 12.0)),
        ],
    )
    append!(gs, A.reaction_gene_association_dnf(model, "ANTI_15378_29101_NhaA")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "ANTI_15378_29101_NhaA")))

    model.reactions["ANTI_15378_29101_NhaB"] = Reaction(
        name = "Antiporter Na/H (NhaB)",
        stoichiometry = Dict(
            "15378" => -3.0, # H+
            "29101_p" => -2.0, # Na+
            "29101" => 2.0, # Na+
            "15378_p" => 3.0, # H+
        ),
        objective_coefficient = 0.0,
        lower_bound = 0,
        upper_bound = 0,
        gene_association = [
            X.Isozyme(; gene_product_stoichiometry = Dict("PN96_RS03465" .=> 12.0)),
        ],
    )
    append!(gs, A.reaction_gene_association_dnf(model, "ANTI_15378_29101_NhaB")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "ANTI_15378_29101_NhaB")))

    model.reactions["ANTI_15378_29101_NhaC"] = Reaction( # stoichiometry is unknown
        name = "Antiporter Na/H (NhaC)",
        stoichiometry = Dict(
            "15378" => -1.0, # H+
            "29101_p" => -1.0, # Na+
            "29101" => 1.0, # Na+
            "15378_p" => 1.0, # H+
        ),
        objective_coefficient = 0.0,
        lower_bound = -1000,
        upper_bound = 1000,
        gene_association = [
            X.Isozyme(; gene_product_stoichiometry = Dict("PN96_RS03215" .=> 12.0)),
            X.Isozyme(; gene_product_stoichiometry = Dict("PN96_RS03265" .=> 12.0)),
            X.Isozyme(; gene_product_stoichiometry = Dict("PN96_RS04845" .=> 12.0)),
            X.Isozyme(; gene_product_stoichiometry = Dict("PN96_RS05055" .=> 12.0)),
            X.Isozyme(; gene_product_stoichiometry = Dict("PN96_RS07920" .=> 12.0)),
            X.Isozyme(; gene_product_stoichiometry = Dict("PN96_RS10500" .=> 12.0)),
            X.Isozyme(; gene_product_stoichiometry = Dict("PN96_RS19810" .=> 12.0)),
            X.Isozyme(; gene_product_stoichiometry = Dict("PN96_RS18650" .=> 12.0)),
            X.Isozyme(; gene_product_stoichiometry = Dict("PN96_RS19835" .=> 12.0)),
            X.Isozyme(; gene_product_stoichiometry = Dict("PN96_RS16655" .=> 12.0)),
        ],
    )
    append!(gs, A.reaction_gene_association_dnf(model, "ANTI_15378_29101_NhaC")...)
    append!(ms, keys(A.reaction_stoichiometry(model, "ANTI_15378_29101_NhaC")))

    add_genes!(model, gs)
    # no need to add metabolites, because they should all already be in the model
    @assert all(in.(ms, Ref(A.metabolites(model))))

end

