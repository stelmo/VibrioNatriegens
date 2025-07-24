
"""
$(TYPEDSIGNATURES)

Build the model.
"""
function build_model()
    #! NB MUST BE THE REFERENCE RHEA REACTION
    df = DataFrame(
        CSV.File(joinpath(pkgdir(@__MODULE__), "data", "model", "metabolic_reactions.csv")),
    )
    get_reactions(unique(df.RHEA_ID)) # per-cache everything

    heteros = @rsubset(df, !ismissing(:Isozyme))
    @select!(heteros, :RHEA_ID, :Protein, :Stoichiometry, :Isozyme)

    homos = @rsubset(df, ismissing(:Isozyme))
    @select!(homos, :RHEA_ID, :Protein, :Stoichiometry)

    ghomos = groupby(homos, [:RHEA_ID, :Protein])
    gheteros = groupby(heteros, [:RHEA_ID, :Isozyme])

    # Build model

    model = Model()

    extend_model!(model, ghomos)
    extend_model!(model, gheteros)
    gapfill!(model)
    curate!(model) # needs to happen here, as metabolites get added that are necessary later

    # # these are all manually added reactions
    add_exchanges!(model)
    add_periplasm_transporters!(model)
    add_membrane_transporters!(model)
    add_electron_transport_chain!(model)
    add_salt_transducers!(model)
    add_atpm!(model)
    add_biomass!(model)

    add_gene_annotations!(model)
    add_metabolite_annotations!(model)
    add_reaction_annotations!(model)

    set_default_exchanges!(model)
    name_reactions!(model)
    name_genes!(model)
    shortname_reactions!(model)

    haskey(model.reactions, "PERM_29101") && delete!(model.reactions, "PERM_29101") # house keeping for Na+
    model
end
