
using RheaReactions, CSV, DataFrames, DataFramesMeta
using AbstractFBCModels
import AbstractFBCModels as A
import COBREXA as X
using DocStringExtensions
import SparseArrays as S
using JSON
using FASTX


gene_seqs = Dict{String, String}()

FASTAReader(open(joinpath("data", "genome", "proteome.fasta"))) do reader
    for record in reader
        pid = string(first(split(string(last(split(description(record), "protein_id="))), "]")))
        length(pid) <= 20 && begin
            gene_seqs[pid] = sequence(record) 
        end
    end
end

seq_df = DataFrame(ProteinAccession=collect(keys(gene_seqs)), Sequence=collect(values(gene_seqs)))

gene_df = DataFrame(
    CSV.File(
        joinpath(
            "data",
            "annotations",
            "ncbi",
            "refseq_annotations.tsv",
        ),
    ),
)
@rename!(gene_df, :ProteinAccession = $"Protein accession", :GeneID = $"Gene ID")
@select!(gene_df, :Name, :ProteinAccession, :Symbol, :Begin, :End, :Chromosome, :Orientation, :Accession, :GeneID)

# secondary place to lookup gene symbol
gene_df2 = DataFrame(
    CSV.File(
        joinpath(
            "data",
            "annotations",
            "eggnog",
            "out.emapper.annotations",
        ),
    ),
)
@select!(gene_df2, :query, :Preferred_name)
@rename!(gene_df2, :ProteinAccession = :query)
leftjoin!(gene_df, gene_df2, on=:ProteinAccession, matchmissing=:notequal)
def_symbol(s1, s2) = begin
    ismissing(s1) || return s1
    ismissing(s2) && return missing
    s2 == "-" && return missing
    return s2
end
@transform!(gene_df, :Symbol = def_symbol.(:Symbol, :Preferred_name))
@subset!(gene_df, .!ismissing.(:ProteinAccession))

molar_masses = DataFrame(
    CSV.File(
        joinpath(
            "data",
            "genome",
            "vnat_molar_masses.csv",
        ),
    )
)
@rename!(molar_masses, :ProteinAccession = :Protein)

leftjoin!(gene_df, molar_masses, on=:ProteinAccession)

leftjoin!(gene_df, seq_df, on=:ProteinAccession)
CSV.write("gene_annotations.csv", gene_df)

# metabolites

names_df = DataFrame(
    CSV.File(
        joinpath(
            "data",
            "chebi",
            "names.tsv",
        ),
    ),
)
@subset!(names_df, :LANGUAGE .== "en")
names_dict = Dict(
    "CHEBI:"*string(first(gdf.COMPOUND_ID)) =>
        unique(gdf.NAME)
     for gdf in groupby(names_df, :COMPOUND_ID)
)

inchi_df = DataFrame(
    CSV.File(
        joinpath(
            "data",
            "chebi",
            "chebiId_inchi.tsv",
        ),
    ),
)
inchi_dict = Dict(zip("CHEBI:".*string.(inchi_df.CHEBI_ID), inchi_df.InChI))

chemical_df = DataFrame(
    CSV.File(
        joinpath(
            "data",
            "chebi",
            "chemical_data.tsv",
        ),
    ),
)
@subset!(chemical_df, :SOURCE .== "ChEBI", :TYPE .== "MASS")
chemical_dict = Dict(
    "CHEBI:"*string(first(gdf.COMPOUND_ID)) => parse(Float64, first(gdf.CHEMICAL_DATA))
     for gdf in groupby(chemical_df, :COMPOUND_ID)
)

accession_df = DataFrame(
    CSV.File(
        joinpath(
            "data",
            "chebi",
            "database_accession.tsv",
        ),
    ),
)
kegg_compound = Dict(
    "CHEBI:"*string(first(gdf.COMPOUND_ID)) => gdf.ACCESSION_NUMBER
    for gdf in groupby(@subset(accession_df, :TYPE .== "KEGG COMPOUND accession"), :COMPOUND_ID)
)

df = DataFrame(
    MetaboliteID = Union{Missing, String}[], 
    KeGG = Union{Missing, Vector{String}}[], 
    MolarMass = Union{Missing, Float64}[], 
    InChI = Union{Missing, String}[], 
    Names = Union{Missing, Vector{String}}[],
)
for k in keys(inchi_dict)
    push!(
        df, 
        (k, get(kegg_compound, k, missing), get(chemical_dict, k, missing), get(inchi_dict, k, missing), get(names_dict, k, missing))
    )
end
df
dropmissing!(df)
CSV.write(joinpath("data", "annotations", "reduced", "metabolite_annotations.csv"),df)
