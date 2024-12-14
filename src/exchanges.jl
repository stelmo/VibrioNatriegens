"""
$(TYPEDSIGNATURES)

Add exchange reactions for basic metabolites
"""
function add_exchanges!(model)

    substrates = [
        "C00221" # beta_glucose
        "C00059" # so4
        "C00007" # o2
        "C00014" # nh4
        "C00697" # n2
        "C00244" # nitrate
        "C00121" # ribose
        "C00672" # 2-Deoxy-D-ribose 1-phosphate
        "C00257" # D-Gluconic acid
        "C01801" # Deoxyribose
        "C00259" # L arabinose
        "C01721" # L-Fuculose
        "C00310" # D-Xylulose 
        "C00216" # d-arabinose
        "C00508" # L-ribulose
        "C00312" # L-xylulose
        "C00507" # L-Rhamnose
        "C01019" # L-fucose
        "C00095" # D-fructose
        "C00124" # D-galactose
        "C00243" # lactose
        "C05796" # galactan
        "C00492" # raffinose
        "C01613" # stachyose
        "C00208" # maltose
        "C00329" # glucosamine
        "C00140" # n acetyl glucosamine
        "C01674" # Chitobiose
        "C00275" # mannose 6 phosphate
        "C00546" # methylglyoxal
        "C05993" # acetyl adenylate
        "C00988" # Phosphoglycolic acid
        "C00263" # homoserine
        "C00154" # Palmitoyl-coa
        "C00116" # glycerol
    ]

    products = [
        "C00267" # alpha glucose
        "C00031" # D-Glucose
        "C00620" # alpha ribose 1 phosphate
        "C00022" # pyruvate
        "C00033" # acetate
        "C00058" # formate
        "C00186" # (S)-Lactate
        "C00469" # ethanol
        "C00122" # fumarate
        "C00042" # succinate
        "C00024" # Acetyl-CoA
        "C00119" # PRPP    
        "C00026" # 2-Oxoglutarate
        "C00027" # h2o2
        "C00167" # UDP-glucuronate
        "C01904" # D-Arabitol
        "C00111" # glycerone-P
        "C00266" # glycoaldehyde
        "C00424" # L-lactaldehyde
        "C01050" # UDP-N-acetylmuramate 
        "C00498" # ADP-glucose (to starch)
        "C00636" # Mannose 1 phosphate
        "C00256" # R lactate
        "C00160" # glycolate
        "C00099" # beta alanine
        "C00168" # Hydroxypyruvate
        "C06010" # (S)-2-Acetolactate
        "C00877" # crotonyl-coa
        "C00288" # hco3
        "C00064" # L-glutamine
        "C00025" # glutamate
        "C00283" # hydrogen sulfide
        "C00097" # L-cysteine
        "C02084" # Tetrathionate
        "C00087" # sulfur
        "C06423" # Octanoic acid
        "C01571" # Decanoic acid
        "C02679" # Dodecanoic acid
        "C06424" # Tetradecanoic acid
        "C00154" # Palmitoyl-CoA
        "C00249" # Hexadecanoic acid
        "C01530" # Octadecanoic acid
        "C00641" # 1,2-Diacyl-sn-glycerol
    ]

    bidirs = [
        "C00282" # h2
        "C00075" # UTP
        "C00013" # diphosphate
        "C00004" # NADH
        "C00003" # NAD+
        "C00005" # NADPH
        "C00006" # NADP+
        "C00011" # CO2
        "C00001" # water
        "C00020" # AMP
        "C00002" # atp
        "C00009" # pi
        "C00008" # adp
        "C00080" # h
        "C00010" # coa
        "C15602" # Quinone
        "C15603" # Hydroquinone
        "C00016" # FAD
        "C01352" # FADH2
        "C00138" # Reduced ferredoxin
        "C00139" # Oxidized ferredoxin
        "C00037" # glycine
        "C00065" # serine
        "C00163" # propanoate
        "C00100" # propanoyl coa
        "C00342" # thioredoxin
        "C00343" # thioredoxin disulfide
        "C00125" # ferricytochrome c
        "C00126" # ferrocytochrome c
        "C00924" # Ferrocytochrome 
        "C00923" # Ferricytochrome 
        "C00101" # Tetrahydrofolate 
        "C00143" # 5,10-Methylenetetrahydrofolate
        "C00234" # 10-Formyltetrahydrofolate
        "C00332" # acetoacetyl-coa
        "C01127" # 4-Hydroxy-2-oxoglutarate
        "C05984" # 2-Hydroxybutanoic acid
        "C00222" # 3-Oxopropanoate
        "C00041" # L-alanine
        "C05359" # electron
        "C00054" # Adenosine 3',5'-bisphosphate (PAP)
        "C00229" # ACP
    ]

    for mid in substrates
        mid in A.metabolites(model) || continue
        nm = A.metabolite_name(model, mid)
        model.reactions["EX_"*mid] = Reaction(;
            name = "Exchange $nm",
            stoichiometry = Dict(mid => -1),
            lower_bound = -10.0,
            upper_bound = 0,
        )
    end

    for mid in products
        mid in A.metabolites(model) || continue
        nm = A.metabolite_name(model, mid)
        model.reactions["EX_"*mid] = Reaction(;
            name = "Exchange $nm",
            stoichiometry = Dict(mid => -1),
            lower_bound = 0,
            upper_bound = 1000,
        )
    end

    for mid in bidirs
        mid in A.metabolites(model) || continue
        nm = A.metabolite_name(model, mid)
        model.reactions["EX_"*mid] = Reaction(;
            name = "Exchange $nm",
            stoichiometry = Dict(mid => -1),
            lower_bound = -1000,
            upper_bound = 1000,
        )
    end

end
