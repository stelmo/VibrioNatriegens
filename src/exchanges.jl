"""
$(TYPEDSIGNATURES)

Add exchange reactions for basic metabolites
"""
function add_exchanges!(model)

    substrates = [
        "C00221" # beta_glucose
        # "C00267" # alpha glucose
        # "C00031", # D-Glucose
    ]

    products = [
        "C00001" # water
        "C00007" # o2
        "C00009" # pi
        "C00011" # CO2
        "C00013" # Diphosphate
        "C00014" # nh4
        "C00022" # pyruvate
        "C00024" # Acetyl-CoA
        "C00033" # acetate
        "C00058" # formate
        "C00059" # so4
        "C00080" # H+
        "C00186" # (S)-Lactate
        "C00469" # ethanol
    ]

    bidirs = [
        "C00004" # NADH
        "C00003" # NAD+
        "C00005" # NADPH
        "C00006" # NADP+
        "C00068" # Thiamin diphosphate
        # "C00020" # AMP
    ]

    for mid in substrates
        mid in A.metabolites(model) || continue
        model.reactions["EX_"*mid] = Reaction(;
            name = "Exchange $mid",
            stoichiometry = Dict(mid => -1),
            lower_bound = -10.0,
            upper_bound = 0,
        )
    end

    for mid in products
        mid in A.metabolites(model) || continue
        model.reactions["EX_"*mid] = Reaction(;
            name = "Exchange $mid",
            stoichiometry = Dict(mid => -1),
            lower_bound = 0,
            upper_bound = 1000,
        )
    end
    
    for mid in bidirs
        mid in A.metabolites(model) || continue
        model.reactions["EX_"*mid] = Reaction(;
            name = "Exchange $mid",
            stoichiometry = Dict(mid => -1),
            lower_bound = -1000,
            upper_bound = 1000,
        )
    end
    
end
