
function add_biomass!(model)

    rows = CSV.File(
        joinpath(pkgdir(@__MODULE__), "data", "model", "biomass.csv");
        types = [String, Float64],
    )
    biomass = Dict(rows.Metabolite .=> rows.Coefficient)

    # required atp for growth hydrolysis equation
    atp_req = 71.06
    # atp_req = 0.0

    biomass["30616"] = biomass["30616"] - atp_req # atp
    biomass["15377"] = -atp_req # water
    biomass["43474"] = atp_req # pi
    biomass["456216"] = atp_req # adp
    biomass["15378"] = atp_req # h+

    model.reactions["BIOMASS"] = Reaction(
        name = "Biomass reaction",
        stoichiometry = biomass,
        objective_coefficient = 1.0,
        lower_bound = 0,
        upper_bound = 1000,
        annotations = Dict("SBO" => ["SBO_0000629"], "acronym" => ["BIOMASS", "biomass"]),
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
        lower_bound = 83.22,
        # lower_bound = 0.0,
        upper_bound = 1000,
        annotations = Dict("SBO" => ["SBO_0000630"], "acronym" => ["ATPM", "atpm"]),
    )
end
