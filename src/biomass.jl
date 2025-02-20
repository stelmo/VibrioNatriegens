"""
$(TYPEDSIGNATURES)

Add a biomass function
"""
function add_biomass!(model)
    model.reactions["biomass"] = Reaction(
        name = "Biomass reaction",
        stoichiometry = Dict(
            # "C00186" => -1, # lactate
            # "C00469" => -1, # ethanol
            # "C00033" => -1, # acetate
            # "C00058" => -1, # formate
            # "C00022" => -1, # pyruvate
            # "C00122" => -1, # fumarate
            # "C00042" => -1, # succinate

            # "C00044" => -1, # GTP
            # "C00002" => -1, # ATP
            # "C00063" => -1, # CTP
            # "C00075" => -1, # UTP

            # "C00286" => -1, # dGTP
            # "C00131" => -1, # dATP
            # "C00458" => -1, # dCTP
            # "C00459" => -1, # dTTP

            # "C00188" => -1, # ectoine

            # "C00041" => -1, # alanine
            # "C00049" => -1, # aspartate
            # "C00152" => -1, # asparagine
            # "C00025" => -1, # glutamate
            # "C00064" => -1, # glutamine
            # "C00188" => -1, # threonine
            # "C00065" => -1, # serine
            # "C00097" => -1, # cysteine
            # "C00037" => -1, # glycine
            # "C00073" => -1, # methionine
            # "C00407" => -1, # isoleucine
            # "C00183" => -1, # valine
            # "C00123" => -1, # leucine
            # "C00047" => -1, # lysine
            # "C00062" => -1, # arginine
            # "C00135" => -1, # histidine
            # "C00082" => -1, # tyrosine
            # "C00079" => -1, # phenylalanine
            # "C00078" => -1, # tryptophan
            # "C00148" => -1, # proline

            # "C06423" => -1, # octanoic acid
            # "C01571" => -1, # decanoic acid
            # "C02679" => -1, # dodecanoic acid
            # "C06424" => -1, # tetradecanoic acid
            # "C00249" => -1, # hexadecanoic acid
            # "C01530" => -1, # octadecanoic acid
        ),
        objective_coefficient = 1.0,
        lower_bound = 0,
    )
end

function add_atpm!(model)
    model.reactions["ATPM"] = Reaction(
        name = "ATP Maintenance reaction",
        stoichiometry = Dict(
            "CHEBI:30616" => -1, # atp
            "CHEBI:15377" => -1, # water
            "CHEBI:43474" => 1, # pi
            "CHEBI:456216" => 1, # adp
            "CHEBI:15378" => 1,  # h+
        ),
        objective_coefficient = 1.0,
        lower_bound = 0,
        upper_bound = 1000,
    )
end
