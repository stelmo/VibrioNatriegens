"""
$(TYPEDSIGNATURES)

Add a biomass function
"""
function add_biomass!(model)
    st = JSON.parsefile(joinpath(pkgdir(@__MODULE__), "data", "model", "biomass.json"))
    model.reactions["biomass"] = Reaction(
        name = "Biomass reaction",
        stoichiometry = Dict(keys(st) .=> float.(values(st))),
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
            "CHEBI:30616" => -1, # atp
            "CHEBI:15377" => -1, # water
            "CHEBI:43474" => 1, # pi
            "CHEBI:456216" => 1, # adp
            "CHEBI:15378" => 1,  # h+
        ),
        objective_coefficient = 0.0,
        lower_bound = 8.0,
        upper_bound = 1000,
        annotations = Dict("SBO" => ["SBO_0000630"]),
    )
end
