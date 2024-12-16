using DataFrames, DataFramesMeta, XLSX, CSV

df = DataFrame(
    XLSX.readtable(
        joinpath("data", "curation", "curated", "base_reactions.xlsx"),
        "metabolism",
    ),
)

find_transport(x, y, port) = begin
    ismissing(x) && return false
    occursin(port, x) || occursin(port, y)
end

anno = DataFrame(CSV.File(joinpath("data", "curation", "all_annotations.csv")))
transporters = @rsubset(anno, find_transport(:KEGG_Description, :RefSeq_Description, "Na+"))
@select!(
    transporters,
    :EN_Reaction,
    :Protein,
    :KEGG_Description,
    :HAMAP_Subunit,
    :Uniprot_Subunit,
    :EN_Symbol,
    :RefSeq_Description,
    :EN_TC,
)
CSV.write("salt.csv", @orderby(transporters, :EN_TC, :KEGG_Description))

rxns = filter(
    x -> x ∉ unique(df.Reaction),
    VibrioNatriegens.get_kegg_reactions_in_pathway("map02010"),
)
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
