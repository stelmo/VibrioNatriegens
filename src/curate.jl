
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

rename(model, rid, nm) = begin
    model.reactions[rid].name = nm
end

function curate!(model)

    bi_dir(model, "10231")
    bi_dir(model, "18796")
    for_dir(model, "18133")
    for_dir(model, "30734")
    for_dir(model, "31210")
    for_dir(model, "30946")
    for_dir(model, "28357")
    for_dir(model, "31042")
    for_dir(model, "28361")
    for_dir(model, "30750")

    # add these metabolites manually

    id = "CHEBI:29101" # Na+
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

    # modify rhea reactions
#     model.reactions["38215"].stoichiometry[]

#     delete!(model.metabolites)
#     38215,"-1.0*NAD(+) + -1.0*D-glucose 6-phosphate <-> 1.0*H(+) + 1.0*NADH + 1.0*6-phospho-D-glucono-1,5-lactone"

# 15841,"-1.0*D-glucose 6-phosphate + -1.0*NADP(+) <-> 1.0*H(+) + 1.0*NADPH + 1.0*6-phospho-D-glucono-1,5-lactone"

# 17825,-1.0*ATP + -1.0*D-glucose -> 1.0*H(+) + 1.0*D-glucose 6-phosphate + 1.0*ADP


    
end

