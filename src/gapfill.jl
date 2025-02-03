
function gapfill!(model)

    rhea_rids = [
        54528 # aldehydo-D-ribose 5-phosphate <-> D-ribose 5-phosphate
        21387
        41935
        50947
        17972
        30122
        43439
        31254
        19048 # pyruvate to citramalate, c5 branched dibasic metabolism
        30134
        30118
        30126
        19740
        32630
        31462
        45139
        22564
        15652 # Spontaneous
        28237 # Spontaneous
        22491 # formaldehyde detox

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

