
function gapfill!(model)
    #! NB MUST BE THE REFERENCE RHEA REACTION
    rhea_rids = [
        54528 # aldehydo-D-ribose 5-phosphate <-> D-ribose 5-phosphate
        19045 # pyruvate to citramalate, c5 branched dibasic metabolism
        15649 # Spontaneous
        28234 # Spontaneous
        22488 # formaldehyde detox
        21368 # urate to 5-hydroxyisourate
        17029 # s-allantoin to allantoate
        63204 # release acp octadecanoate
        41932 # release acp hexadecanoate
        30123 # release acp tetradecanoate
        30119 # release acp dodecanoate
        30115 # release acp decanoate
        30131 # release acp octanoate
        19045 # pyruvate to citramalate in valine, leucine, isoleucine biosynthesis
        25078 # Spontaneous
        18609 # arginine biosynthesis
        19533 # arginine & proline metabolism
        16953 # arginine & proline metabolism
        19737 # arginine & proline metabolism
        28234 # Spontaneous
        13389 # tryptophan metabolism
        21384 # phenylalanine metabolism
        12921 # cofactors, make 5-aminolevulinate, actually hemL should be active, but requires tRNA reactions
        54644 # d-amino acid metabolism, connect S)-1-pyrroline-5-carboxylate to make trans-4-hydroxy-L-proline
        14277 # rhamnose metabolism degrade lactaldehyde
    ]
    get_reactions(rhea_rids)

    n = length(rhea_rids)
    df = DataFrame(
        RHEA_ID = rhea_rids,
        Protein = fill(nothing, n),
        Stoichiometry = fill(1),
        Isozyme = fill("A", n),
    )
    dfs = groupby(df, :RHEA_ID)

    extend_model!(model, dfs)
end

