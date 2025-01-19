
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
    for_dir(model, "18136")
    for_dir(model, "30734")
    for_dir(model, "31210")
    for_dir(model, "30946")
    for_dir(model, "28357")
    for_dir(model, "31042")
    for_dir(model, "28361")
    for_dir(model, "30750")

    # add these metabolites 

    id = "CHEBI:29101" # Na+
    model.metabolites[id] = Metabolite(
        name = "Na+",
        formula = Dict("Na" => 1),
        compartment = "Cytosol",
        charge = 1,
    )

    
end
