
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

rhea_rxn_dir(rxn, qrt) = begin 
    idx = first(indexin([rxn], qrt))
    isnothing(idx) && error("Reaction not found...")
    idx == 1 && return (-1000, 1000)
    idx == 2 && return (0, 1000)
    idx == 3 && return (-1000, 0)
    idx == 4 && return (-1000, 1000)
end

function curate!(model)

    # add these metabolites manually

    id = "CHEBI:29101" # Na+ this is needed since the salt reactions only get added in the transporter section
    model.metabolites[id] = Metabolite(
        name = "Na+",
        formula = Dict("Na" => 1),
        compartment = "Cytosol",
        charge = 1,
    )

    id = "RHEA-COMP:14399" # Fe(III)-[cytochrome c]
    model.metabolites[id] = Metabolite(
        name = "Fe(III)-[cytochrome c]",
        formula = Dict("Fe" => 1),
        compartment = "Cytosol",
        charge = 3,
    )

    id = "RHEA-COMP:10350" # Fe(II)-[cytochrome c]
    model.metabolites[id] = Metabolite(
        name = "Fe(II)-[cytochrome c]",
        formula = Dict("Fe" => 1),
        compartment = "Cytosol",
        charge = 2,
    )

    id = "CHEBI:17992" # sucrose
    model.metabolites[id] = Metabolite(
        name = "Sucrose",
        formula = Dict("C" => 12,"H" => 22,"O" => 11,),
        compartment = "Cytosol",
        charge = 0,
    )

    # adjust the formula of [thioredoxin]-dithiol: S1C3N1H5O1 -> S2C6N2H10O2
    model.metabolites["CHEBI:29950"].formula = Dict(
        "S" => 2,
        "C" => 6,
        "N" => 2,
        "H" => 10,
        "O" => 2,
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

    # change directions to match what is found in biocyc - manual thermodynamics leaves much to be desired
    biocyc = DataFrame(CSV.File(joinpath("data", "annotations", "rhea", "biocyc_rxns.csv")))
    @select!(biocyc, :rheaDir, :metacyc)
    for rid in A.reactions(model)
        qrt = RheaReactions.get_reaction_quartet(parse(Int, rid))
        df = @subset(biocyc, in.(:rheaDir, Ref(qrt)))
        isempty(df) && continue
        lb, ub = rhea_rxn_dir(df[1,1], qrt)
        model.reactions[rid].lower_bound = lb
        model.reactions[rid].upper_bound = ub    
    end

    # change directions manually
    bi_dir(model, "13196") # make GTP, dATP, dGTP work
    for_dir(model, "20309") # h2o2 -> o2 + h2o
    for_dir(model, "19028") # h2o2 producer
    for_dir(model, "30778") # h2o2 producer
    for_dir(model, "15056") # prevent loop in nucleotides through 15056 <-> 31134
    for_dir(model, "31134") # prevent loop in nucleotides through 15056 <-> 31134 
    for_dir(model, "41812") # prevent loop in lipids through 41812 <-> 54868
    for_dir(model, "54868") # prevent loop in lipids through 41812 <-> 54868 
    for_dir(model, "41528") # prevent loop in lipids through 41528 <-> 41848
    for_dir(model, "41848") # prevent loop in lipids through 41528 <-> 41848 
    for_dir(model, "54936") # prevent loop in lipids through 54936 <-> 41864
    for_dir(model, "41864") # prevent loop in lipids through 54936 <-> 41864 
    for_dir(model, "54900") # prevent loop in lipids through 54900 <-> 41912
    for_dir(model, "41912") # prevent loop in lipids through 54900 <-> 41912 
    for_dir(model, "15308") # prevent loop in lipids through 15308 <-> 30070
    for_dir(model, "30070") # prevent loop in lipids through 15308 <-> 30070 
    for_dir(model, "38215") # prevent loop in carbohydrates through 38215 <-> 38215
    for_dir(model, "38215") # prevent loop in carbohydrates through 38215 <-> 38215 
    

end

