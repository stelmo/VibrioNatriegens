using CSV, DataFrames, DataFramesMeta, HTTP, VibrioNatriegens

eggnog = DataFrame(
    CSV.File(
        joinpath("data", "annotations", "eggnog", "out.emapper.annotations");
        missingstring = ["", "-"],
    ),
)
@select!(eggnog, :query, :EC, :KEGG_ko, :KEGG_Reaction, :KEGG_TC, :Preferred_name)
@rename!(
    eggnog,
    :Protein = :query,
    :EN_KO = :KEGG_ko,
    :EN_Reaction = :KEGG_Reaction,
    :EN_TC = :KEGG_TC,
    :EN_EC = :EC,
    :EN_Symbol = :Preferred_name
)
@rtransform!(eggnog, :EN_KO = ismissing(:EN_KO) ? missing : split(:EN_KO, ","))
eggnog = flatten(eggnog, :EN_KO; scalar = Missing)
@rtransform!(eggnog, :EN_KO = ismissing(:EN_KO) ? missing : String(:EN_KO[4:end]))

kegg = DataFrame(
    CSV.File(
        joinpath("data", "annotations", "kegg", "ko.txt"),
        header = ["Protein", "KO", "Description"],
    ),
)
@rename!(kegg, :KEGG_KO = :KO, :KEGG_Description = :Description)
dropmissing!(kegg)
anno = leftjoin(eggnog, kegg, on = :Protein, matchmissing = :notequal)

hamap = DataFrame(CSV.File(joinpath("data", "annotations", "hamap", "hamap_subunits.csv")))
@rename!(hamap, :HAMAP_Subunit = :Subunit)
leftjoin!(anno, hamap, on = :Protein)

refseq =
    DataFrame(CSV.File(joinpath("data", "annotations", "ncbi", "refseq_annotations.tsv")))
@rename!(refseq, :Protein = $"Protein accession")
@select!(refseq, :Protein, :Symbol, :Name)
@rename!(refseq, :RefSeq_Symbol = :Symbol, :RefSeq_Description = :Name)
leftjoin!(anno, refseq, on = :Protein, matchmissing = :notequal)

ecoli =
    DataFrame(CSV.File(joinpath("data", "ecoli", "e_coli_gene_names_subunit_stoich.tsv")))
@rename!(
    ecoli,
    :Uniprot_Symbol = $"Gene Names (primary)",
    :Uniprot_Subunit = $"Subunit structure",
    :Uniprot_Accession = :Entry
)
@select!(ecoli, :Uniprot_Subunit, :Uniprot_Symbol, :Uniprot_Accession)
dropmissing!(ecoli)
@rtransform!(ecoli, :Uniprot_Subunit = :Uniprot_Subunit[10:end])

leftjoin!(anno, ecoli, on = (:EN_Symbol => :Uniprot_Symbol), matchmissing = :notequal)
@rtransform!(
    anno,
    :EN_Reaction = ismissing(:EN_Reaction) ? missing : split(:EN_Reaction, ",")
)
anno = flatten(anno, :EN_Reaction; scalar = Missing)

kr(id) = begin # need to catch here, because some eggnog reactions are invalid
    println(id)
    try
        VibrioNatriegens.get_kegg_reaction(string(id)).string_stoichiometry
    catch err
        missing
    end
end

kegg_rxns = collect(skipmissing(anno.EN_Reaction))
kegg_strs = kr.(kegg_rxns)

rid_str_lu = Dict(kegg_rxns .=> kegg_strs)

@rtransform!(anno, :KEGG_Reaction_Definition = get(rid_str_lu, :EN_Reaction, missing))
anno.KEGG_Reaction_Definition

CSV.write(joinpath("data", "curation", "all_annotations.csv"), anno)
