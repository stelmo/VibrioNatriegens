"""
$(TYPEDSIGNATURES)

Add a biomass function
"""
function add_biomass!(model)
    model.reactions["biomass"] = Reaction(
        name = "Biomass reaction",
        stoichiometry = Dict(
            "C00186" => -1, # lactate
            # "C00469" => -1, # ethanol
            # "C00033" => -1, # acetate
            # "C00058" => -1, # formate
            # "C00022" => -1, # pyruvate
            # "C00122" => -1, # fumarate
            # "C00042" => -1, # succinate
            # "C00119" => -1, # PRPP
            # "C00044" => -1, # GTP
            # "C00002" => -1, # ATP
            # "C00063" => -1, # CTP
            # "C00075" => -1, # UTP
            # "C00286" => -1, # dGTP
            # "C00131" => -1, # dATP
            # "C00458" => -1, # dCTP
            # "C00459" => -1, # dTTP
            # "C00041" => -1, # alanine
            # "C00049" => -1, # aspartate
            # "C00152" => -1, # asparagine
            # "C00025" => -1, # glutamate
            # "C00064" => -1, # glutamine
            # "C00188" => -1, # threonine
            # "C00188" => -1, # ectoine
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
            # "C00378" => -1, # thiamine
            # "C00068" => -1, # thiamine diphosphate
            # "C00016" => -1, # FAD
            # "C00255" => -1, # riboflavin
            # # "C00250" => -1, # pyrodoxine (vit B6)
            # "C00003" => -1, # NAD
            # "C00006" => -1, # NADP
            # "C00864" => -1, # pantothenate
            # "C00120" => -1, # biotin
            # "C00504" => -1, # folate
            # "C00143" => -1, # 5,10 methylene THF
            # "C00234" => -1, # 10 formyl THF
            # "C15670" => -1, # heme A
            # "C00390" => -1, # ubiquinol
            # "C06423" => -1, # octanoic acid
            # "C01571" => -1, # decanoic acid
            # "C02679" => -1, # dodecanoic acid
            # "C06424" => -1, # tetradecanoic acid
            # "C00249" => -1, # hexadecanoic acid
            # "C01530" => -1, # octadecanoic acid
            # # "C0" => -1, # 
            # "C0" => -1, # 
            # "C0" => -1, # 
            # "C0" => -1, # 
            # "C0" => -1, # 
            # "C0" => -1, # 
            # "C0" => -1, # 
        ),
        objective_coefficient = 1.0,
        lower_bound = 0,
    )
end

function add_atpm!(model)
    model.reactions["ATPM"] = Reaction(
        name = "ATP Maintenance reaction",
        stoichiometry = Dict(
            "C00002" => -1, # atp
            "C00001" => -1, # water
            "C00009" => 1, # pi
            "C00008" => 1, # adp
            "C00080" => 1,  # h+
        ),
        objective_coefficient = 0.0,
        lower_bound = 0,
        upper_bound = 1000,
    )
end
