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
    rename(
        model,
        "R03270",
        "pyruvate:thiamin diphosphate acetaldehydetransferase (decarboxylating) [step 2]",
    )
    rename(
        model,
        "R00014",
        "pyruvate:thiamin diphosphate acetaldehydetransferase (decarboxylating) [step 1]",
    )
    bi_dir(model, "R01325")
    bi_dir(model, "R01899")
    bi_dir(model, "R01057")
    rename(model, "R10997", "2-oxoisovalerate dehydrogenase [step 2]")
    rename(model, "R10996", "2-oxoisovalerate dehydrogenase [step 1]")
    rename(model, "R10150", "tetrathionate reductase")
    for_dir(model, "R00529") # dg error
    rename(
        model,
        "R07762",
        "hexadecanoyl-[acyl-carrier protein]:malonyl-[acyl-carrier-protein] C-acyltransferase (decarboxylating)",
    )
    rename(
        model,
        "R07763",
        "oxostearoyl-[acyl-carrier protein]:malonyl-[acyl-carrier-protein] C-acyltransferase (decarboxylating)",
    )
    rename(
        model,
        "R07764",
        "hydroxyoctadecanoyl-[acyl-carrier protein]:malonyl-[acyl-carrier-protein] C-acyltransferase (decarboxylating)",
    )
    for_dir(model, "R01175")
    for_dir(model, "R04751")
    for_dir(model, "R03777")
    for_dir(model, "R04754")
    for_dir(model, "R03857")
    for_dir(model, "R03990")
    for_dir(model, "R01175")
    for_dir(model, "R01279")
    rename(model, "R07390", "Cardiolipin synthase")
    rename(model, "R03270", "2-oxoglutarate dehydrogenase complex [step 2]")
    rename(model, "R03634", "Stachyose hydrolase")

    rename(model, "R02508", "Cystathionine transferase")
    rename(model, "R02670", "Cinnavalininate catalase")
    rename(model, "R03314", "L-Glutamate hydrolase")
    rename(model, "R03634", "Stachyose hydrolase")
    rename(model, "R04326", "5'-Phosphoribosylglycinamide transferase")
    rename(model, "R11264", "2-methylaconitate isomerase")
    rename(model, "R10151", "L-cysteine S-thiosulfotransferase")
    rename(model, "R09840", "phenylacetyl-CoA thioesterase")
    rename(model, "R07603", "Oxo-acid dehydrogenase complex")
    rename(model, "R07396", "4-Methylthio-2-oxobutanoic acid transferase")
    rename(model, "R04546", "Propanoyl-CoA oxidation")
    rename(
        model,
        "R04594",
        "N1-(5-Phospho-alpha-D-ribosyl)-5,6-dimethylbenzimidazole hydrolase",
    )
    rename(model, "R05221", "Adenosyl cobinamide transferase")
    rename(model, "R05223", "Adenosine-GDP-cobinamide transferase")
    rename(model, "R05225", "adenosylcobyric acid synthase")
    rename(model, "R05237", "trans-3-Chloroallyl aldehyde dehydrogenase (NAD+)")
    rename(model, "R05238", "cis-3-Chloroallyl aldehyde dehydrogenase (NAD+)")
    rename(model, "R05290", "benzoate/toluate 1,2-dioxygenase")
    rename(model, "R05428", "toluate dioxygenase")
    rename(model, "R05506", "benzoyl acetyl-CoA thiolase")
    rename(model, "R05549", "alpha-galactosidase")
    rename(model, "R05586", "3-Oxopimeloyl-CoA acetyl-CoA acyltransferase")
    rename(model, "R05595", "Crotonoyl-CoA hydratase")
    rename(model, "R05621", "benzoate/toluate 1,2-dioxygenase")
    rename(model, "R05665", "benzoate/toluate 1,2-dioxygenase")
    rename(
        model,
        "R06558",
        "adenosylcobinamide kinase / adenosylcobinamide-phosphate guanylyltransferase",
    )
    rename(model, "R07316", "Carbamate lyase")
    rename(model, "R07399", "(R)-citramalate synthase")
    rename(model, "R07406", "4-hydroxythreonine-4-phosphate dehydrogenase")
    rename(model, "R07415", "gamma-glutamylputrescine oxidase")
    rename(model, "R07599", "2-oxoisovalerate dehydrogenase")
    rename(model, "R07600", "2-oxoisovalerate dehydrogenase")
    rename(model, "R07601", "2-oxoisovalerate dehydrogenase")
    rename(model, "R07602", "2-oxoisovalerate dehydrogenase")
    rename(model, "R07604", "2-oxoisovalerate dehydrogenase")
    rename(model, "R07891", "OPC6-coa acetyl-CoA acyltransferase")
    rename(model, "R07895", "OPC4-coa acetyl-CoA acyltransferase")
    rename(model, "R07899", "Isojasmonic acid transferase")

end
