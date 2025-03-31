"""
$(TYPEDSIGNATURES)

Add a biomass function
"""
function add_biomass!(model)
    biomass = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "model", "biomass.json"))

    # required atp for growth hydrolysis equation
    atp_req = 110.5380339572664
    # atp_req = 0.0
    
    biomass["30616"] = biomass["30616"] - atp_req # atp
    biomass["15377"] = -atp_req # water
    biomass["43474"] = atp_req # pi
    biomass["456216"] = atp_req # adp
    biomass["15378"] = atp_req # h+


    model.reactions["biomass"] = Reaction(
        name = "Biomass reaction",
        stoichiometry = Dict(keys(biomass) .=> float.(values(biomass))),
        objective_coefficient = 1.0,
        lower_bound = 0,
        upper_bound = 1000,
        annotations = Dict("SBO" => ["SBO_0000629"]),
    )
end

function add_atpm!(model)
    model.reactions["ATPM"] = Reaction(
        name = "ATP Maintenance reaction",
        stoichiometry = Dict(
            "30616" => -1, # atp
            "15377" => -1, # water
            "43474" => 1, # pi
            "456216" => 1, # adp
            "15378" => 1,  # h+
        ),
        objective_coefficient = 0.0,
        lower_bound = 71.43699059080485,
        # lower_bound = 0.0,
        upper_bound = 1000,
        annotations = Dict("SBO" => ["SBO_0000630"]),
    )
end
