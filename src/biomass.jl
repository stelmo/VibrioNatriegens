"""
$(TYPEDSIGNATURES)

Add a biomass function
"""
function add_biomass!(model)
    model.reactions["biomass"] = Reaction(
        name = "Biomass reaction",
        stoichiometry = Dict(
            "glycogen" => -1,
            
            "CHEBI:37565" => -1, # GTP
            "CHEBI:37563" => -1, # CTP
            "CHEBI:46398" => -1, # UTP

            "CHEBI:30616" => -1, # atp
            "CHEBI:15377" => -1, # water
            "CHEBI:43474" => 1, # pi
            "CHEBI:456216" => 1, # adp
            "CHEBI:15378" => 1,  # h+

            "CHEBI:61429" => -1, # dGTP
            "CHEBI:61404" => -1, # dATP
            "CHEBI:61481" => -1, # dCTP
            "CHEBI:37568" => -1, # dTTP

            "CHEBI:58515" => -1, # ectoine

            "CHEBI:32551" => -1, # lysine
            "CHEBI:57844" => -1, # methionine
            "CHEBI:35235" => -1, # cysteine
            "CHEBI:57972" => -1, # alanine
            "CHEBI:29991" => -1, # aspartate
            "CHEBI:58048" => -1, # asparagine
            "CHEBI:29985" => -1, # glutamate
            "CHEBI:58359" => -1, # glutamine
            "CHEBI:57926" => -1, # threonine
            "CHEBI:33384" => -1, # serine
            "CHEBI:57305" => -1, # glycine
            "CHEBI:58045" => -1, # isoleucine
            "CHEBI:57762" => -1, # valine
            "CHEBI:57427" => -1, # leucine
            "CHEBI:32682" => -1, # arginine
            "CHEBI:57595" => -1, # histidine
            "CHEBI:58315" => -1, # tyrosine
            "CHEBI:58095" => -1, # phenylalanine
            "CHEBI:57912" => -1, # tryptophan
            "CHEBI:60039" => -1, # proline

            "CHEBI:25646" => -1, # octanoic acid
            "CHEBI:27689" => -1, # decanoic acid
            "CHEBI:18262" => -1, # dodecanoic acid
            "CHEBI:30807" => -1, # tetradecanoic acid
            "CHEBI:7896" => -1, # hexadecanoic acid
            "CHEBI:25629" => -1, # octadecanoic acid
        ),
        objective_coefficient = 1.0,
        lower_bound = 0,
        upper_bound = 1000,
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
        objective_coefficient = 0.0,
        lower_bound = 1,
        upper_bound = 1000,
    )
end
