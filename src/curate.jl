
bi_dir(model, rid) = begin
    model.reactions[rid].lower_bound = -1000.0
    model.reactions[rid].upper_bound = 1000.0
end

for_dir(model, rid) = begin
    model.reactions[rid].lower_bound = 0.0
    model.reactions[rid].upper_bound = 1000.0
end

rev_dir(model, rid) = begin
    model.reactions[rid].lower_bound = -1000.0
    model.reactions[rid].upper_bound = 0.0
end

function curate!(model)

    # bi_dir(model, "10231")
    # bi_dir(model, "18796")
    # for_dir(model, "18133")
    # for_dir(model, "30734")
    # for_dir(model, "31210")
    # for_dir(model, "30946")
    # for_dir(model, "28357")
    # for_dir(model, "31042")
    # for_dir(model, "28361")
    # for_dir(model, "30750")

    # add these metabolites manually

    id = "CHEBI:29101" # Na+ this is needed since the salt reactions only get added in the transporter section
    model.metabolites[id] = Metabolite(
        name = "Na+",
        formula = Dict("Na" => 1),
        compartment = "Cytosol",
        charge = 1,
    )

    id = "CHEBI:15983" # Ferrocytochrome
    model.metabolites[id] = Metabolite(
        name = "Ferrocytochrome",
        formula = Dict("X" => 1),
        compartment = "Cytosol",
        charge = 0,
    )

    id = "CHEBI:15719" # Ferricytochrome
    model.metabolites[id] = Metabolite(
        name = "Ferricytochrome",
        formula = Dict("X" => 1),
        compartment = "Cytosol",
        charge = 0,
    )

    # modify rhea reactions to exchange D- with the beta-D isomer
    delete!(model.metabolites, "CHEBI:4167") # D-glucose
    delete!(model.metabolites, "CHEBI:61548") # D-glucopyranose 6-phosphate

    model.reactions["38215"].stoichiometry["CHEBI:58247"] = model.reactions["38215"].stoichiometry["CHEBI:61548"] # D-glucopyranose 6-phosphate -> 	β-D-glucose 6-phosphate
    delete!(model.reactions["38215"].stoichiometry, "CHEBI:61548")

    model.reactions["15841"].stoichiometry["CHEBI:58247"] = model.reactions["15841"].stoichiometry["CHEBI:61548"] # D-glucopyranose 6-phosphate -> 	β-D-glucose 6-phosphate
    delete!(model.reactions["15841"].stoichiometry, "CHEBI:61548")

    model.reactions["17825"].stoichiometry["CHEBI:15903"] = model.reactions["17825"].stoichiometry["CHEBI:4167"] # D-glucose -> beta-D-glucose
    delete!(model.reactions["17825"].stoichiometry, "CHEBI:4167")

    model.reactions["17825"].stoichiometry["CHEBI:58247"] = model.reactions["17825"].stoichiometry["CHEBI:61548"] # D-glucopyranose 6-phosphate -> β-D-glucose 6-phosphate
    delete!(model.reactions["17825"].stoichiometry, "CHEBI:61548")

    model.reactions["76650"].stoichiometry["CHEBI:71044"] = model.reactions["76650"].stoichiometry["CHEBI:195329"] # (E) 2,3-didehydroadipoyl-CoA -> 2,3-didehydroadipoyl-CoA
    delete!(model.reactions["76650"].stoichiometry, "CHEBI:195329")
    delete!(model.metabolites, "CHEBI:195329")


end

