
# "C02084" # Tetrathionate
# "C00087" # sulfur
# "C06423" # Octanoic acid
# "C01571" # Decanoic acid
# "C02679" # Dodecanoic acid
# "C06424" # Tetradecanoic acid
# "C00154" # Palmitoyl-CoA
# "C00249" # Hexadecanoic acid
# "C01530" # Octadecanoic acid
# "C00641" # 1,2-Diacyl-sn-glycerol
# "C00075" # UTP
# "C00013" # diphosphate
# "C00004" # NADH
# "C00003" # NAD+
# "C00005" # NADPH
# "C00006" # NADP+
# "C00672" # 2-Deoxy-D-ribose 1-phosphate
# "C01801" # Deoxyribose
# "C00020" # AMP
# "C00002" # atp
# "C00008" # adp
# "C00275" # mannose 6 phosphate
# "C00546" # methylglyoxal
# "C05993" # acetyl adenylate
# "C00988" # Phosphoglycolic acid
# "C00263" # homoserine
# "C00154" # Palmitoyl-coa
# "C00620" # alpha ribose 1 phosphate
# "C00024" # Acetyl-CoA
# "C00119" # PRPP    
# "C00026" # 2-Oxoglutarate
# "C00027" # h2o2
# "C00167" # UDP-glucuronate
# "C00111" # glycerone-P
# "C00266" # glycoaldehyde
# "C00424" # L-lactaldehyde
# "C01050" # UDP-N-acetylmuramate 
# "C00498" # ADP-glucose (to starch)
# "C00636" # Mannose 1 phosphate
# "C00160" # glycolate
# "C00099" # beta alanine
# "C00168" # Hydroxypyruvate
# "C06010" # (S)-2-Acetolactate
# "C00877" # crotonyl-coa
# "C05359" # electron
# "C15602" # Quinone
# "C15603" # Hydroquinone
# "C00016" # FAD
# "C01352" # FADH2
# "C00004" # NADH
# "C00003" # NAD+
# "C00138" # Reduced ferredoxin
# "C00139" # Oxidized ferredoxin
# "C00163" # propanoate
# "C00100" # propanoyl coa
# "C00342" # thioredoxin
# "C00343" # thioredoxin disulfide
# "C00125" # ferricytochrome c
# "C00126" # ferrocytochrome c
# "C00924" # Ferrocytochrome 
# "C00923" # Ferricytochrome 
# "C00101" # Tetrahydrofolate 
# "C00143" # 5,10-Methylenetetrahydrofolate
# "C00234" # 10-Formyltetrahydrofolate
# "C00332" # acetoacetyl-coa
# "C01127" # 4-Hydroxy-2-oxoglutarate
# "C05984" # 2-Hydroxybutanoic acid
# "C00222" # 3-Oxopropanoate
# "C00054" # Adenosine 3',5'-bisphosphate (PAP)
# "C00229" # ACP
# "C00010" # coa
# "C00697" # n2
# "C00244" # nitrate
# "C00121" # ribose
# "C00257" # D-Gluconic acid
# "C00259" # L-arabinose
# "C01721" # L-Fuculose
# "C00310" # D-Xylulose 
# "C00216" # D-Arabinose
# "C00508" # L-ribulose
# "C00312" # L-xylulose
# "C00507" # L-Rhamnose
# "C01019" # L-fucose
# "C00095" # D-fructose
# "C00124" # D-galactose
# "C00243" # Lactose
# "C05796" # Galactan
# "C00492" # Raffinose
# "C01613" # Stachyose
# "C00208" # Maltose
# "C00329" # Glucosamine
# "C00140" # N-Acetyl-Glucosamine
# "C01674" # Chitobiose
# "C00116" # Glycerol
# "C00267" # Alpha-glucose
# "C00031" # D-Glucose
# "C00022" # Pyruvate
# "C00033" # Acetate
# "C00058" # Formate
# "C00186" # (S)-Lactate
# "C00469" # Ethanol
# "C00122" # Fumarate
# "C00042" # Succinate
# "C01904" # D-Arabitol
# "C00256" # (R)-lactate
# "C00288" # hco3
# "C00064" # L-glutamine
# "C00025" # glutamate
# "C00283" # hydrogen sulfide
# "C00097" # L-cysteine
# "C00282" # h2
# "C00011" # CO2
# "C00001" # water
# "C00080" # h
# "C00037" # glycine
# "C00065" # serine
# "C00041" # L-alanine

"""
$(TYPEDSIGNATURES)

Add exchange reactions for basic metabolites
"""
function add_exchanges!(model)

    df = DataFrame(
        XLSX.readtable(
            joinpath("data", "curation", "curated", "base_reactions.xlsx"),
            "exchanges",
        ),
    )

    all_exchange_metabolites = df.KeGG

    substrates = [ # allowed to be imported
        "C00221" # beta_glucose
        "C00059" # so4
        "C00007" # o2
        "C00014" # nh4
        "C00009" # pi
    ]

    for mid in all_exchange_metabolites
        mid in A.metabolites(model) || (@warn "Metabolite $mid not in model!"; continue)

        if mid == "C00221" # default carbon source
            lb, ub = (-22.0, 0.0)
        elseif mid in substrates
            lb, ub = (-1000.0, 0.0)
        else # product
            lb, ub = (0.0, 1000.0)
        end
        nm = A.metabolite_name(model, mid)

        model.metabolites[mid*"_e"] = deepcopy(model.metabolites[mid])

        model.metabolites[mid*"_e"].compartment = "Extracellular"

        model.reactions["EX_"*mid] = Reaction(;
            name = "Exchange $nm",
            stoichiometry = Dict(mid * "_e" => -1),
            lower_bound = lb,
            upper_bound = ub,
        )
    end

end

function add_periplasm_transporters!(model)

    # all extracellular move with diffusion into periplasm

    df = DataFrame(
        XLSX.readtable(
            joinpath("data", "curation", "curated", "base_reactions.xlsx"),
            "exchanges",
        ),
    )

    # add salt here, should not have an exchange 
    nm = A.metabolite_name(model, "C01330")
    model.metabolites["C01330_p"] = deepcopy(model.metabolites["C01330"])
    model.metabolites["C01330_p"].compartment = "Periplasm"

    all_exchange_metabolites = df.KeGG
    for mid in all_exchange_metabolites # Na+
        mid in A.metabolites(model) || (@warn "Metabolite $mid not in model!"; continue)

        nm = A.metabolite_name(model, mid)

        model.metabolites[mid*"_p"] = deepcopy(model.metabolites[mid])

        model.metabolites[mid*"_p"].compartment = "Periplasm"

        model.reactions["DF_"*mid] = Reaction(;
            name = "Diffusion $nm",
            stoichiometry = Dict(mid * "_e" => -1, mid * "_p" => 1),
            lower_bound = -1000.0,
            upper_bound = 1000.0,
            gene_association = [
                X.Isozyme(; gene_product_stoichiometry = Dict(["Missing"] .=> [1.0])),
            ],
        )
    end

end

function add_membrane_transporters!(model)
    amino_acids =
        String.(
            DataFrame(
                XLSX.readtable(
                    joinpath("data", "curation", "curated", "base_reactions.xlsx"),
                    "aminoacids",
                ),
            ).KeGG
        )

    df = DataFrame(
        XLSX.readtable(
            joinpath("data", "curation", "curated", "base_reactions.xlsx"),
            "transporters",
        ),
    )
    @select!(df, :Type, :Metabolite, :KeggID, :Protein, :Stoichiometry, :Subunit)
    gs = String[]
    ms = String[]

    # these amino acids ABS are a special cases
    abc_aa = @subset(df, ismissing.(:KeggID))
    general_aa = @subset(abc_aa, :Metabolite .== "AminoAcid")
    branched_aa = @subset(abc_aa, :Metabolite .== "BranchedAminoAcid") # L, I, V

    for g in groupby(general_aa, :Subunit)
        for mid in amino_acids
            push!(ms, mid)
            rid = "ABC_$mid"
            iso = string.(g.Protein)
            ss = parse.(Float64, string.(g.Stoichiometry))
            append!(gs, iso)
            add_abc!(model, rid, mid, iso, ss)
        end
    end

    for g in groupby(branched_aa, :Subunit)
        for mid in ["C00183", "C00407", "C00123"]
            rid = "ABC_$mid"
            push!(ms, mid)
            iso = string.(g.Protein)
            ss = parse.(Float64, string.(g.Stoichiometry))
            append!(gs, iso)
            add_abc!(model, rid, mid, iso, ss)
        end
    end

    # other abc transporters (conventionally added to source)
    dropmissing!(df)
    abcs = @subset(df, :Type .== "ABC")
    for g in groupby(abcs, [:KeggID, :Subunit])
        mid = first(g.KeggID)
        if mid in A.metabolites(model)
            push!(ms, mid)
            rid = "ABC_$mid"
            iso = string.(g.Protein)
            ss = parse.(Float64, string.(g.Stoichiometry))
            append!(gs, iso)
            add_abc!(model, rid, mid, iso, ss)
        else
            @warn "$mid not in model (ABC)"
        end
    end

    # PTS transporters
    pts = @subset(df, :Type .== "PTS")
    for g in groupby(pts, [:KeggID, :Subunit])
        mid = first(g.KeggID)
        if mid in A.metabolites(model)
            push!(ms, mid)
            rid = "PTS_$mid"
            iso = string.(g.Protein)
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_pts!(model, rid, mid, iso, ss)
        else
            @warn "$mid not in model (PTS)"
        end
    end

    # symport
    symport = @subset(df, :Type .== "Symport")
    for g in groupby(symport, [:KeggID, :Subunit])
        mid1, mid2 = sort(split(first(g.KeggID), "/")) # to make rid unique
        if mid1 in A.metabolites(model) && mid2 in A.metabolites(model)
            push!(ms, mid1)
            push!(ms, mid2)
            rid = "Symport_$(mid1)_$mid2"
            iso = string.(g.Protein)
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_symport!(model, rid, mid1, mid2, iso, ss)
        else
            @warn "$mid1 or $mid2 not in model (symport)"
        end
    end

    # antiport
    antiport = @subset(df, :Type .== "Antiport")
    for g in groupby(antiport, [:KeggID, :Subunit])
        mid1, mid2 = sort(split(first(g.KeggID), "/")) # to make rid unique
        if mid1 in A.metabolites(model) && mid2 in A.metabolites(model)
            rid = "Antiport_$(mid1)_$mid2"
            push!(ms, mid1)
            push!(ms, mid2)
            iso = string.(g.Protein)
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_antiport!(model, rid, mid1, mid2, iso, ss)
        else
            @warn "$mid1 or $mid2 not in model (antiport)"
        end
    end

    # permease (the default as well)
    permease = @subset(df, :Type .== "Permease")
    for g in groupby(permease, [:KeggID, :Subunit])
        mid = first(g.KeggID)
        if mid in A.metabolites(model)
            rid = "Permease_$mid"
            iso = string.(g.Protein)
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_permease!(model, rid, mid, iso, ss)
        else
            @warn "$mid not in model (permease)"
        end
    end

    # add default permeases
    all_exchange_metabolites =
        string.(
            DataFrame(
                XLSX.readtable(
                    joinpath("data", "curation", "curated", "base_reactions.xlsx"),
                    "exchanges",
                ),
            ).KeGG
        )
    missing_transporters = setdiff(all_exchange_metabolites, unique(ms))
    for mid in missing_transporters
        if mid in A.metabolites(model)
            rid = "Permease_$mid"
            add_permease!(model, rid, mid, ["Missing"], [1.0])
        else
            @warn "$mid not in model (missing permease)"
        end
    end
end

function add_abc!(model, rid, mid, iso, ss)
    isoz = X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Transport $mid ABC",
            stoichiometry = Dict(
                "C00002" => -1, # atp
                "C00001" => -1, # water
                "C00009" => 1, # pi
                "C00008" => 1, # adp
                "C00080" => 1,  # h+ # uncharged reaction
                mid * "_p" => -1.0,
                mid => 1.0, # cytosol
            ),
            objective_coefficient = 0.0,
            lower_bound = -1000, # can secrete at a cost
            upper_bound = 1000,
            gene_association = [isoz],
        )
    end
end


function add_pts!(model, rid, mid, iso, ss)
    isoz = X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Transport $mid PTS",
            stoichiometry = Dict(
                "C00074" => -1.0, # pep
                "C00022" => 1.0, # pyr
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

function add_symport!(model, rid, mid1, mid2, iso, ss)
    isoz = X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Transport $mid1, $mid2 Symport",
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

function add_antiport!(model, rid, mid1, mid2, iso, ss)
    isoz = X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Transport $mid1, $mid2 Antiport",
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

function add_permease!(model, rid, mid, iso, ss)
    isoz = X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Transport $mid Permease",
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
            "C00138" => -2.0, # Reduced ferredoxin
            "C00004" => -1.0,  # nadh
            "C00006" => -2.0, # NADP+
            "C00080" => -1.0, # H+
            "C00139" => 2.0, # Oxidized ferredoxin
            "C00003" => 1.0, # NAD
            "C00005" => 2.0, # NADPH
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
            "C00002" => 1, # atp
            "C00001" => 1, # water
            "C00009" => -1, # pi
            "C00008" => -1, # adp
            "C00080" => 3.0,  # h+
            "C00080_p" => -4.0,  # h+
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
            "C15602" => 1.0, # Quinone
            "C15603" => -1.0, # Hydroquinone
            "C00126" => 2.0, # Ferrocytochrome c
            "C00125" => -2.0, # Ferricytochrome c
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
            "C00007" => -1, # o2
            "C00126" => 4.0, # Ferrocytochrome c
            "C00080" => -8.0, # H+
            "C00125" => 4.0, # Ferricytochrome c
            "C00001" => 2.0, # h2o
            "C00080_p" => 4.0, # H+
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
            "C00007" => -1, # o2
            "C00126" => 4.0, # Ferrocytochrome c
            "C00080" => -8.0, # H+
            "C00125" => 4.0, # Ferricytochrome c
            "C00001" => 2.0, # h2o
            "C00080_p" => 4.0, # H+
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
            "C00001" => 1.0, # h2o
            "C00080" => -2.0, # H+
            "C00080_p" => 2.0, # H+
            "C00007" => -0.5, # o2
            "C15602" => 1.0, # Quinone
            "C15603" => -1.0, # Hydroquinone
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
            "C00001" => -1.0, # h2o
            "C00080" => -4.0, # H+
            "C00080_p" => 4.0, # H+
            "C00007" => -0.5, # o2
            "C15602" => 1.0, # Quinone
            "C15603" => -1.0, # Hydroquinone
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
            "C00003" => -1.0, # NAD
            "C00005" => -1.0, # NADPH
            "C00080" => -1.0, # H+
            "C00004" => 1.0, # nadh
            "C00006" => 1.0, # NADP+
            "C00080_p" => 1.0, # H+
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
            "C00002" => 1, # atp
            "C00001" => 1, # water
            "C00009" => -1, # pi
            "C00008" => -1, # adp
            "C00080" => -1.0,  # h+
            "C01330_p" => -4.0, # Na+
            "C01330" => 4.0, # Na+
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
            "C00036" => -1.0, # oxaloacetate
            "C01330" => -2.0, # Na+
            "C00080" => -1.0, # H+
            "C00022" => 1.0, # pyruvate
            "C01330_p" => 2.0, # Na+
            "C00011" => 1.0, # CO2
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
            "C01330" => 1.0, # Na+ (n?)
            "C00003" => 1.0, # nad
            "C15603" => 1.0, # Hydroquinone
            "C15602" => -1.0, # Quinone
            "C01330_p" => -1.0, # Na+ (n?)
            "C00004" => -1.0, # nadh
            "C00080" => -1.0, # h+
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
            "C00138" => -2, # Reduced ferredoxin
            "C01330" => -1, # Na+
            "C00003" => -1, # NAD
            "C00080" => -1, # H+
            "C00139" => 2, # Oxidized ferredoxin
            "C01330_p" => 1.0, # Na+
            "C00004" => 1.0,  # nadh
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



