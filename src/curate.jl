
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

function set_default_exchanges()


    substrates = [ # allowed to be imported
    "CHEBI:4167" # glucose
    "CHEBI:16189" # so4
    "CHEBI:15379" # o2
    "CHEBI:28938" # nh4(+)
    "CHEBI:43474" # pi
    "CHEBI:29101" # Na+
]

bidirs = [
    "CHEBI:15377" # H2O
]

for mid in String.(df.CHEBI)

    if mid == "CHEBI:4167" # default carbon source
        lb, ub = (-22.0, 0.0)
    elseif mid in substrates
        lb, ub = (-1000.0, 0.0)
    elseif min in bidirs
        lb, ub = (-1000.0, 1000.0)
    else
        lb, ub = (0.0, 1000.0)
    end
    
end
