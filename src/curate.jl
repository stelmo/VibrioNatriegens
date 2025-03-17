
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

rhea_rxn_dir(rid, consensus) = begin
    idx = consensus - parse(Int, rid)
    idx == 0 && return (-1000, 1000)
    idx == 1 && return (0, 1000)
    idx == 2 && return (-1000, 0)
    idx == 3 && return (-1000, 1000)
end

function curate!(model)

    # delete reactions added to load specific metabolites in gapfill
    delete!(model.reactions, "19289") # sucrose
    delete!(model.reactions, "27814") # Na+ 

    # glycogen
    id = "glycogen"
    model.metabolites[id] = Metabolite(
        name = "Glycogen",
        formula = Dict("C" => 6.0, "H" => 10.0, "O" => 5.0),
        compartment = "Cytosol",
        charge = 0,
    )

    # modify rhea reactions to exchange D- with the beta-D isomer
    delete!(model.metabolites, "4167") # D-glucose
    delete!(model.metabolites, "61548") # D-glucopyranose 6-phosphate
    delete!(model.metabolites, "75989") # alpha-D-glucosamine 6-phosphate
    delete!(model.metabolites, "195329") # (E) 2,3-didehydroadipoyl-CoA 

    model.reactions["12172"].stoichiometry["58725"] =
        model.reactions["12172"].stoichiometry["75989"] # alpha-D-glucosamine 6-phosphate -> D-glucosamine 6-phosphate
    delete!(model.reactions["12172"].stoichiometry, "75989")

    model.reactions["38215"].stoichiometry["58247"] =
        model.reactions["38215"].stoichiometry["61548"] # D-glucopyranose 6-phosphate -> 	β-D-glucose 6-phosphate
    delete!(model.reactions["38215"].stoichiometry, "61548")

    model.reactions["15841"].stoichiometry["58247"] =
        model.reactions["15841"].stoichiometry["61548"] # D-glucopyranose 6-phosphate -> 	β-D-glucose 6-phosphate
    delete!(model.reactions["15841"].stoichiometry, "61548")

    model.reactions["17825"].stoichiometry["15903"] =
        model.reactions["17825"].stoichiometry["4167"] # D-glucose -> beta-D-glucose
    delete!(model.reactions["17825"].stoichiometry, "4167")

    model.reactions["17825"].stoichiometry["58247"] =
        model.reactions["17825"].stoichiometry["61548"] # D-glucopyranose 6-phosphate -> β-D-glucose 6-phosphate
    delete!(model.reactions["17825"].stoichiometry, "61548")

    model.reactions["76647"].stoichiometry["71044"] =
        model.reactions["76647"].stoichiometry["195329"] # (E) 2,3-didehydroadipoyl-CoA -> 2,3-didehydroadipoyl-CoA
    delete!(model.reactions["76647"].stoichiometry, "195329")

    model.reactions["33791"].stoichiometry["37721"] =
        model.reactions["33791"].stoichiometry["28645"] # β-D-Fructose (28645) -> D-Fructose (37721), we want it to be D-Fructose, RHEA:33791 bc sucrose hydrolase
    delete!(model.reactions["33791"].stoichiometry, "28645")

    # change directions to match what is found in biocyc - manual thermodynamics leaves much to be desired

    metacyc = DataFrame(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "rhea", "biocyc_rxns.csv"),
            drop = [3],
        ),
    )
    @rename!(metacyc, :metacyc = :rheaDir)

    ecocyc = DataFrame(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "rhea", "ecocyc_rxns.csv"),
            drop = [3],
        ),
    )
    @rename!(ecocyc, :ecocyc = :rheaDir)

    j = outerjoin(metacyc, ecocyc, on = :rhea)
    @rtransform!(
        j,
        :consensus = ismissing(:ecocyc) ? :metacyc : :ecocyc,
        :rhea = string(:rhea)
    )
    @subset!(j, in.(:rhea, Ref(A.reactions(model))))

    for (rid, consensus) in zip(j.rhea, j.consensus)
        lb, ub = rhea_rxn_dir(rid, consensus)
        model.reactions[rid].lower_bound = lb
        model.reactions[rid].upper_bound = ub
    end

    # change directions manually
    bi_dir(model, "13193") # make GTP, dATP, dGTP work
    for_dir(model, "20309") # h2o2 -> o2 + h2o
    for_dir(model, "30775") # h2o2 producer
    for_dir(model, "15053") # prevent loop in nucleotides through 15056 <-> 31134
    for_dir(model, "31131") # prevent loop in nucleotides through 15056 <-> 31134 
    for_dir(model, "41812") # prevent loop in lipids through 41812 <-> 54868
    for_dir(model, "54868") # prevent loop in lipids through 41812 <-> 54868 
    for_dir(model, "41528") # prevent loop in lipids through 41528 <-> 41848
    for_dir(model, "41848") # prevent loop in lipids through 41528 <-> 41848 
    for_dir(model, "54936") # prevent loop in lipids through 54936 <-> 41864
    for_dir(model, "41864") # prevent loop in lipids through 54936 <-> 41864 
    for_dir(model, "54900") # prevent loop in lipids through 54900 <-> 41912
    for_dir(model, "41912") # prevent loop in lipids through 54900 <-> 41912 
    for_dir(model, "15305") # prevent loop in lipids through 15308 <-> 30070
    for_dir(model, "30067") # prevent loop in lipids through 15308 <-> 30070 
    for_dir(model, "38215") # prevent loop in carbohydrates through 38215 <-> 38215 
    rev_dir(model, "11612") # prevent loop + metacyc GLUTAMATESYN-RXN
    for_dir(model, "11692") # https://www.uniprot.org/uniprotkb/P27306/entry
    for_dir(model, "27758") # https://biocyc.org/reaction?orgid=ECOLI&id=GCVMULTI-RXN
    for_dir(model, "24793") # histidine biosynthesis blocked, https://biocyc.org/reaction?orgid=META&id=GLUTAMIDOTRANS-RXN (incorrect ref in rhea wrt direction)    
    rev_dir(model, "17737") # F, Y, W biosynthesis blocked, https://biocyc.org/reaction?orgid=META&id=SHIKIMATE-5-DEHYDROGENASE-RXN (rhea has both directions listed)
    rev_dir(model, "25290") # prevent loop in ethanol production 
    rev_dir(model, "22852") # serine biosyn natural dir + loop prevention
    rev_dir(model, "19125") # serine biosyn natural dir + loop prevention
    for_dir(model, "18793") # loop prev
    for_dir(model, "29391") # loop prev
    for_dir(model, "23352") # glutamate synthesis, loop prev 
    for_dir(model, "32263") # 4 aminobut degrad, loop prev 
    rev_dir(model, "30799") # loop prev
    rev_dir(model, "16197") # loop prev
    for_dir(model, "13213") # loop prev
    for_dir(model, "13217") # loop prev
    rev_dir(model, "15869") # d alanine syn only + loop prev
    rev_dir(model, "15133") # glut biosyn
    
    # add custom reactions (needs to happen after direction setting)

    model.reactions["glycogen_synthase"] = Reaction(
        name = "Glycogen synthase (ADPGlc)",
        stoichiometry = Dict(
            "57498" => -1.0, # ADP-alpha-D-glucose
            "456216" => 1.0, # adp
            "glycogen" => 1.0, # 
            "15378" => 1.0, # h+
        ),
        objective_coefficient = 0.0,
        lower_bound = 0,
        upper_bound = 1000,
        gene_association = [
            X.Isozyme(; gene_product_stoichiometry = Dict("WP_020332873.1" .=> 1.0)),
        ],
    )

    model.reactions["glycogen_phosphorylase"] = Reaction(
        name = "Glycogen phosphorylase",
        stoichiometry = Dict(
            "glycogen" => -1.0,
            "58601" => 1.0, # alpha-D-glucose 1-phosphate
            "43474" => -1.0, # phosphate 
        ),
        objective_coefficient = 0.0,
        lower_bound = 0,
        upper_bound = 1000,
        gene_association = [
            X.Isozyme(; gene_product_stoichiometry = Dict("WP_020333475.1" .=> 1.0)),
        ],
    )
    add_genes!(model, ["WP_020333475.1", "WP_020332873.1"])

    model.reactions["polyphosphate_kinase2"] = Reaction(
        name = "ppk2",
        stoichiometry = Dict(
            "456216" => 1.0, # adp
            "30616" => -1.0, # atp
            "33019" => -1.0, # diphospate
            "18036" => 1.0, # triphosphate
            "15378" => 1.0, # h+
        ),
        objective_coefficient = 0.0,
        lower_bound = 0,
        upper_bound = 1000,
        gene_association = [
            X.Isozyme(; gene_product_stoichiometry = Dict("WP_020333630.1" .=> 4.0)),
        ],
    )
    model.reactions["polyphosphate_kinase1"] = Reaction(
        name = "ppk1",
        stoichiometry = Dict(
            "456216" => 1.0, # adp
            "30616" => -1.0, # atp
            "43474" => -1.0, # phosphate 
            "33019" => 1.0, # diphospate
        ),
        objective_coefficient = 0.0,
        lower_bound = 0,
        upper_bound = 1000,
        gene_association = [
            X.Isozyme(; gene_product_stoichiometry = Dict("WP_020333630.1" .=> 4.0)),
        ],
    )
    add_genes!(model, ["WP_020333630.1"])

end

