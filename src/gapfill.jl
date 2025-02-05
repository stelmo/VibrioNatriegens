
function gapfill!(model)

    rhea_rids = [
        54528 # aldehydo-D-ribose 5-phosphate <-> D-ribose 5-phosphate
        19048 # pyruvate to citramalate, c5 branched dibasic metabolism
        15652 # Spontaneous
        28237 # Spontaneous
        22491 # formaldehyde detox
        21371 # urate to 5-hydroxyisourate
        17032 # s-allantoin to allantoate
        # 30118
        # 30122
        # 30126
        # 30134
        # 41935
        21387
        50947
        17972
        43439
        31254
        19740
        32630
        31462
        45139
        # 42703 look at this one, gene yes, should be gone?
    ]


    n = length(rhea_rids)
    df = DataFrame(
        RHEA_ID = rhea_rids,
        Protein = fill("Missing", n),
        Stoichiometry = fill(1, ),
        Subunit = fill("A", n),
        DeltaG = fill(nothing, n),
        RevIndex = fill(nothing, n),
    )
    dfs = groupby(df, :RHEA_ID)

    VibrioNatriegens.extend_model!(model, dfs)
end

