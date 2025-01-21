
"""
$(TYPEDSIGNATURES)

Add exchange reactions for basic metabolites
"""
function add_exchanges!(model)

    df = DataFrame(
        CSV.File(
            joinpath("data", "model", "exchange_metabolites.csv"),
        ),
    )

    # by default, exchanges can only export metabolites
    for mid in String.(df.CHEBI)

        nm = A.metabolite_name(model, mid)

        model.metabolites[mid*"_e"] = deepcopy(model.metabolites[mid])

        model.metabolites[mid*"_e"].compartment = "Extracellular"

        model.reactions["EX_"*mid] = Reaction(;
            name = "Exchange $nm",
            stoichiometry = Dict(mid * "_e" => -1),
            lower_bound = 0,
            upper_bound = 1000,
        )
    end

end

function add_periplasm_transporters!(model)

    # all extracellular move with diffusion into periplasm
    df = DataFrame(
        CSV.File(
            joinpath("data", "model", "exchange_metabolites.csv"),
        ),
    )

    # bidirectional by default
    for mid in String.(df.CHEBI)

        nm = A.metabolite_name(model, mid)

        model.metabolites[mid*"_p"] = deepcopy(model.metabolites[mid])

        model.metabolites[mid*"_p"].compartment = "Periplasm"

        model.reactions["DF_"*mid] = Reaction(;
            name = "Diffusion $nm",
            stoichiometry = Dict(mid * "_e" => -1, mid * "_p" => 1),
            lower_bound = -1000.0,
            upper_bound = 1000.0,
        )
    end

end

function add_membrane_transporters!(model)
   
    df = DataFrame(
        CSV.File(
            joinpath("data", "model", "transporters.csv"),
        ),
    )

    gs = String[]
    ms = String[]

    # abc transporters
    abcs = @subset(df, :Type .== "ABC")
    for g in groupby(abcs, [:CHEBI, :Subunit])
        mid = first(g.CHEBI)
        if mid in A.metabolites(model)
            push!(ms, mid)
            iso = string.(g.Protein)
            ss = parse.(Float64, string.(g.Stoichiometry))
            append!(gs, iso)
            add_abc!(model, mid, iso, ss)
        else
            @warn "$mid not in model (ABC)"
        end
    end

    # PTS transporters
    pts = @subset(df, :Type .== "PTS")
    for g in groupby(pts, [:CHEBI, :Subunit])
        mid = first(g.CHEBI)
        if mid in A.metabolites(model)
            push!(ms, mid)
            iso = string.(g.Protein)
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_pts!(model, mid, iso, ss)
        else
            @warn "$mid not in model (PTS)"
        end
    end

    # symport
    symport = @subset(df, :Type .== "Symport")
    for g in groupby(symport, [:CHEBI, :Subunit])
        mid1, mid2 = sort(split(first(g.CHEBI), "/")) # to make rid unique
        if mid1 in A.metabolites(model) && mid2 in A.metabolites(model)
            push!(ms, mid1)
            push!(ms, mid2)
            iso = string.(g.Protein)
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_symport!(model, mid1, mid2, iso, ss)
        else
            @warn "$mid1 or $mid2 not in model (symport)"
        end
    end

    # antiport
    antiport = @subset(df, :Type .== "Antiport")
    for g in groupby(antiport, [:CHEBI, :Subunit])
        mid1, mid2 = sort(split(first(g.CHEBI), "/")) # to make rid unique
        if mid1 in A.metabolites(model) && mid2 in A.metabolites(model)
            push!(ms, mid1)
            push!(ms, mid2)
            iso = string.(g.Protein)
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_antiport!(model, mid1, mid2, iso, ss)
        else
            @warn "$mid1 or $mid2 not in model (antiport)"
        end
    end

    # permease (the default as well)
    permease = @subset(df, :Type .== "Permease")
    for g in groupby(permease, [:CHEBI, :Subunit])
        mid = first(g.CHEBI)
        if mid in A.metabolites(model)
            push!(ms, mid)
            iso = string.(g.Protein)
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_permease!(model, mid, iso, ss)
        else
            @warn "$mid not in model (permease)"
        end
    end

    # add default permeases
    all_exchange_metabolites =  DataFrame(
        CSV.File(
            joinpath("data", "model", "exchange_metabolites.csv"),
        ),
    ).CHEBI
    # note: Pi, Na, and H will not get a permease here, due to them being involved in the other porters
    missing_transporters = setdiff(all_exchange_metabolites, unique(ms))
    for mid in missing_transporters
        if mid in A.metabolites(model)
            add_permease!(model, mid, ["Missing"], [1.0])
        end
    end
end

function add_abc!(model, mid, iso, ss)
    rid = "ABC_$mid"
    isoz = X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Transport $(A.metabolite_name(model, String(mid))) ABC",
            stoichiometry = Dict(
                "CHEBI:30616" => -1, # atp
                "CHEBI:15377" => -1, # water
                "CHEBI:43474" => 1, # pi
                "CHEBI:456216" => 1, # adp
                "CHEBI:15378" => 1,  # h+ 
                mid * "_p" => -1.0,
                mid => 1.0, # cytosol
            ),
            objective_coefficient = 0.0,
            lower_bound = 0,
            upper_bound = 1000,
            gene_association = [isoz],
        )
    end
end

function add_pts!(model, mid, iso, ss)
    rid = "PTS_$mid"
    isoz = X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Transport $(A.metabolite_name(model, String(mid))) PTS",
            stoichiometry = Dict(
                "CHEBI:58702" => -1.0, # pep
                "CHEBI:15361" => 1.0, # pyr
                mid * "_p" => -1.0,
                mid => 1.0, # cytosol
            ),
            objective_coefficient = 0.0,
            lower_bound = 0,
            upper_bound = 1000,
            gene_association = [isoz],
        )
    end
end

function add_symport!(model, mid1, mid2, iso, ss)
    rid = "SYM_$(mid1)_$mid2"
    isoz = X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Symport $(A.metabolite_name(model, String(mid1)))::$(A.metabolite_name(model, String(mid2)))",
            stoichiometry = Dict(
                mid1 * "_p" => -1.0,
                mid1 => 1.0,
                mid2 * "_p" => -1.0,
                mid2 => 1.0,
            ),
            objective_coefficient = 0.0,
            lower_bound = -1000,
            upper_bound = 1000,
            gene_association = [isoz],
        )
    end
end

function add_antiport!(model, mid1, mid2, iso, ss)
    rid = "ANTI_$(mid1)_$mid2"
    isoz = X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Antiport $(A.metabolite_name(model, String(mid1)))::$(A.metabolite_name(model, String(mid2)))",
            stoichiometry = Dict(
                mid1 * "_p" => 1.0,
                mid1 => -1.0,
                mid2 * "_p" => -1.0,
                mid2 => 1.0,
            ),
            objective_coefficient = 0.0,
            lower_bound = -1000,
            upper_bound = 1000,
            gene_association = [isoz],
        )
    end
end

function add_permease!(model, mid, iso, ss)
    rid = "PERM_$mid"
    isoz = X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Permease $(A.metabolite_name(model, String(mid)))",
            stoichiometry = Dict(mid * "_p" => -1.0, mid => 1.0),
            objective_coefficient = 0.0,
            lower_bound = -1000,
            upper_bound = 1000,
            gene_association = [isoz],
        )
    end
end

function add_electron_transport_chain!(model)

    model.reactions["R-nfn"] = Reaction(; # this is speculative
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
        ),
    )

    model.reactions["R-H-ATPsynthase"] = Reaction(;
        name = "F-type H+-transporting ATPase",
        stoichiometry = Dict(
            "CHEBI:30616" => 1, # atp
            "CHEBI:15377" => 1, # water
            "CHEBI:43474" => -1, # pi
            "CHEBI:456216" => -1, # adp
            "CHEBI:15378" => 3.0,  # h+
            "CHEBI:15378_p" => -4.0,  # h+
        ),
        lower_bound = -1000.0,
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

    model.reactions["R-cyt-bc1"] = Reaction(;
        name = "Cytochrome c oxidase",
        stoichiometry = Dict(
            "CHEBI:132124" => 1.0, # a ubiquinone
            "CHEBI:24646" => -1.0, # a ubiquinol
            "CHEBI:15983" => 2.0, # Ferrocytochrome c
            "CHEBI:15719" => -2.0, # Ferricytochrome c
        ),
        lower_bound = -1000.0,
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
        annotations = Dict(
            "KEGG_REACTION" => [
                "Hydroquinone + 2 Ferricytochrome c <=> Quinone + 2 Ferrocytochrome c + 2 H+",
            ],
            "EC" => ["7.1.1.8"],
        ),
    )

    model.reactions["R-cyt-c"] = Reaction(;
        name = "Cytochrome c oxidase",
        stoichiometry = Dict(
            "CHEBI:15379" => -1, # o2
            "CHEBI:15983" => 4.0, # Ferrocytochrome c
            "CHEBI:15378" => -8.0, # H+
            "CHEBI:15719" => 4.0, # Ferricytochrome c
            "CHEBI:15377" => 2.0, # h2o
            "CHEBI:15378_p" => 4.0, # H+
        ),
        lower_bound = -1000.0,
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
        annotations = Dict(
            "KEGG_REACTION" => [
                "Oxygen + 4 Ferrocytochrome c + 8 H+ <=> 4 Ferricytochrome c + 2 H2O + 4 H+",
            ],
            "EC" => ["7.1.1.9"],
        ),
    )

    model.reactions["R-cyt-c-cbb3"] = Reaction(;
        name = "Cytochrome c oxidase, cbb3-type",
        stoichiometry = Dict(
            "CHEBI:15379" => -1, # o2
            "CHEBI:15983" => 4.0, # Ferrocytochrome c
            "CHEBI:15378" => -8.0, # H+
            "CHEBI:15719" => 4.0, # Ferricytochrome c
            "CHEBI:15377" => 2.0, # h2o
            "CHEBI:15378_p" => 4.0, # H+
        ),
        lower_bound = -1000.0,
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
        annotations = Dict(
            "KEGG_REACTION" => [
                "Oxygen + 4 Ferrocytochrome c + 8 H+ <=> 4 Ferricytochrome c + 2 H2O + 4 H+",
            ],
            "EC" => ["7.1.1.9"],
        ),
    )

    model.reactions["R-cyt-bd"] = Reaction(;
        name = "Cytochrome BD-I",
        stoichiometry = Dict(
            "CHEBI:15377" => 1.0, # h2o
            "CHEBI:15378" => -2.0, # H+
            "CHEBI:15378_p" => 2.0, # H+
            "CHEBI:15379" => -0.5, # o2
            "CHEBI:132124" => 1.0, # a ubiquinone
            "CHEBI:24646" => -1.0, # a ubiquinol
        ),
        lower_bound = -1000.0,
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
        annotations = Dict(
            "KEGG_REACTION" =>
                ["2 Ubiquinol + Oxygen + 4 H+ <=> 2 Ubiquinone + 2 H2O + 4 H+"],
            "EC" => ["7.1.1.7"],
        ),
    )

    model.reactions["R-cyt-bo"] = Reaction(;
        name = "Cytochrome oxidase bo3 (ubiquinol: 4 protons) (periplasm)",
        stoichiometry = Dict(
            "CHEBI:15377" => -1.0, # h2o
            "CHEBI:15378" => -4.0, # H+
            "CHEBI:15378_p" => 4.0, # H+
            "CHEBI:15379" => -0.5, # o2
            "CHEBI:132124" => 1.0, # a ubiquinone
            "CHEBI:24646" => -1.0, # a ubiquinol
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
        annotations = Dict(
            "KEGG_REACTION" => ["h2o + 4 h_c + 0.5 o2 + quinol -> quinone + 4 h_p"],
        ),
    )

    model.reactions["R-pnt"] = Reaction(;
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
        ),
    )


end

function add_salt_transducers!(model)

    model.reactions["R-Na-ATPsynthase"] = Reaction(;
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
        lower_bound = -1000.0,
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
        ),
    )

    model.reactions["R-oad"] = Reaction(;
        name = "oxaloacetate decarboxylase (Na(+) extruding)",
        stoichiometry = Dict(
            "CHEBI:16452" => -1.0, # oxaloacetate
            "CHEBI:29101" => -2.0, # Na+
            "CHEBI:15378" => -1.0, # H+
            "CHEBI:15361" => 1.0, # pyruvate
            "CHEBI:29101_p" => 2.0, # Na+
            "CHEBI:16526" => 1.0, # CO2
        ),
        lower_bound = -1000.0,
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

    model.reactions["R-nqr"] = Reaction(;
        name = "Na+-transporting NADH:ubiquinone oxidoreductase",
        stoichiometry = Dict(
            "CHEBI:29101" => 1.0, # Na+ (n?)
            "CHEBI:57540" => 1.0, # nad
            "CHEBI:24646" => 1.0, # a ubiquinol
            "CHEBI:132124" => -1.0, # a ubiquinone
            "CHEBI:29101_p" => -1.0, # Na+ (n?)
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
            "EXPASY" => ["https://enzyme.expasy.org/EC/7.2.1.1"],
        ),
    )

    model.reactions["R-rnf"] = Reaction(;
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

end



