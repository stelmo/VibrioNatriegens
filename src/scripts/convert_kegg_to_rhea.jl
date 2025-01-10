using DataFrames, CSV, DataFramesMeta, XLSX, Serialization

# load curated list of reactions
df = DataFrame(
    XLSX.readtable(
        joinpath("data", "curation", "curated", "base_reactions.xlsx"),
        "metabolism",
    ),
)

meta = DataFrame(CSV.File(joinpath("data", "metanetx", "reac_xref.tsv")))
rhea = DataFrame(CSV.File(joinpath("data", "rhea", "rhea2kegg_reaction.tsv")))

# first get rhea official cross mapping
@rename!(rhea, :Reaction = :ID)
@select!(rhea, :RHEA_ID, :Reaction)
df_rhea = innerjoin(df, rhea, on = :Reaction)

# join the stragglers
df_rems = @rsubset(df, :Reaction ∉ unique(df_rhea.Reaction))
df_metakegg = @rsubset(meta, occursin("kegg.reaction", :source))
@rtransform!(df_metakegg, :Reaction = last(split(:source, ":")))
@select!(df_metakegg, :ID, :Reaction)
df_rems2 = innerjoin(df_rems, df_metakegg, on = :Reaction)
@select!(df_rems2, :Reaction, :ID)
df_metarhea = @rsubset(meta, occursin("rhea:", :source))
@rsubset!(df_metarhea, occursin("<=>", :description))
df_rems3 = innerjoin(df_metarhea, df_rems2, on = :ID)
@rtransform!(df_rems3, :RHEA_ID = parse(Int, last(split(:source, ":"))))
@select!(df_rems3, :Reaction, :RHEA_ID)
@transform!(df_rems3, :Reaction = string.(:Reaction), :RHEA_ID = :RHEA_ID)
append!(df_rems3, DataFrame(CSV.File(joinpath("data", "curation", "curated", "missing_rids.csv"))))
@combine(groupby(df_rems3, :Reaction), :Reaction = first(:Reaction), :RHEA_ID = first(:RHEA_ID)) # unique matchups only
df_rhea2 = innerjoin(df, df_rems3, on = :Reaction)

# final set
df2 = append!(df_rhea, df_rhea2)
gg(x) = ismissing(x) ? missing : string.(x)
@transform!(
    df2,
    :Reaction = string.(:Reaction),
    :Protein = string.(:Protein),
    :Stoichiometry = string.(:Stoichiometry),
    :Subunit = gg.(:Subunit),
    :Notes = gg.(:Notes)
)

rhea_rxn = deserialize("rhea-reaction.js")
gg(x) = haskey(rhea_rxn, string(x)) ? rhea_rxn[string(x)].def : missing
@rtransform!(df2, :Equation = gg(:RHEA_ID))

old_model = DataFrame(CSV.File(joinpath("data", "curation", "curated", "model.csv")))
@select!(old_model, :Reaction, :Stoichiometry)
@rename!(old_model, :KEGG = :Stoichiometry)

leftjoin!(df2, old_model, on=:Reaction)
CSV.write(joinpath("data", "curation", "curated", "rhea_base.csv"), @orderby(df2, :Reaction))




