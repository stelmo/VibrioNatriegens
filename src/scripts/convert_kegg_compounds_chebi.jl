using XLSX, DataFrames, DataFramesMeta, CSV, VibrioNatriegens

c = VibrioNatriegens._get_cache("compounds", "CHEBI")
r = VibrioNatriegens._get_cache("reactions", "RHEA")
r["21308"]

df = DataFrame(
    XLSX.readtable(
        joinpath("data", "curation", "curated", "base_reactions.xlsx"),
        "exchanges",
    ),
)

chebi = DataFrame(CSV.File(joinpath("data", "chebi", "database_accession.tsv")))
@select!(chebi, :COMPOUND_ID, :SOURCE,:ACCESSION_NUMBER)
@rename!(chebi, :KeGG = :ACCESSION_NUMBER)
@rsubset!(chebi, startswith(:KeGG, "C"))

chebi = @combine(
    groupby(chebi, :KeGG),
    :KeGG = first(:KeGG),
    :CHEBI = first(sort(:COMPOUND_ID)),
)

ph = DataFrame(CSV.File(joinpath("data", "chebi", "chebi_pH7_3_mapping.tsv")))
phlu = Dict(ph.CHEBI .=> ph.CHEBI_PH7_3)
_lu = Dict(chebi.KeGG .=> chebi.CHEBI)

lu = Dict(k => get(phlu, v, v) for (k, v) in _lu )

dff = @rtransform(
    df,
    :KeGG = lu[string(:KeGG)],
)
@rtransform!(
    dff,
    :Name = haskey(c, string(:KeGG)) ? c[string(:KeGG)].entry : "Missing"
)



CSV.write("chebi.csv", dfff)


