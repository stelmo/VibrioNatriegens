using DataFrames, DataFramesMeta, XLSX, CSV

df = DataFrame(
    XLSX.readtable(
        joinpath("data", "curation", "curated", "base_reactions.xlsx"),
        "Sheet1",
    ),
)
rxns = filter(
    x -> x ∉ unique(df.Reaction),
    VibrioNatriegens.get_kegg_reactions_in_pathway("map00780"),
)


anno = DataFrame(CSV.File(joinpath("data", "curation", "all_annotations.csv")))
@rsubset!(anno, !ismissing(:KEGG_Reaction_Definition))
@rsubset!(anno, !ismissing(:EN_KO), !ismissing(:EN_Reaction))
@rename!(anno, :Reaction = :EN_Reaction)

@select!(
    anno,
    :Reaction,
    :Protein,
    :EN_EC,
    :EN_KO,
    :KEGG_KO,
    :HAMAP_Subunit,
    :Uniprot_Subunit,
    :EN_Symbol,
    :KEGG_Description,
    :RefSeq_Description,
    :KEGG_Reaction_Definition,
)

@rsubset!(anno, :Reaction in rxns)
CSV.write("pathway_temp.csv", @orderby(anno, :Reaction))
