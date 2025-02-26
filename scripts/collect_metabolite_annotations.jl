using CSV, DataFrames, DataFramesMeta
using VibrioNatriegens, JSON
import AbstractFBCModels as A

model = VibrioNatriegens.build_model()
mids =
    last.(
        split.(
            filter(
                x -> !occursin("_", x),
                filter(startswith("CHEBI"), A.metabolites(model)),
            ),
            ":",
        )
    )
mids = parse.(Int, mids)

chemical_data = DataFrame(CSV.File(joinpath("data", "chebi", "chemical_data.tsv")))
ms = unique(chemical_data.COMPOUND_ID)
cd = DataFrame(
    ChebiID = String[],
    Charge = Union{String,Missing}[],
    Mass = Union{Missing,String}[],
    Formula = Union{Missing,String}[],
)
for df in groupby(chemical_data, :COMPOUND_ID)
    d = Dict(zip(df.TYPE, df.CHEMICAL_DATA))
    push!(
        cd,
        (
            "CHEBI:" * string(first(df.COMPOUND_ID)),
            get(d, "CHARGE", missing),
            get(d, "MASS", missing),
            get(d, "FORMULA", missing),
        ),
    )
end


inchi = DataFrame(CSV.File(joinpath("data", "chebi", "chebiId_inchi.tsv")))

accession_data = DataFrame(CSV.File(joinpath("data", "chebi", "database_accession.tsv")))
groupby(accession_data, :COMPOUND_ID)

