
"""
$(TYPEDSIGNATURES)

Build the model.
"""
function build_model()

    df = DataFrame(
        CSV.File(joinpath(pkgdir(@__MODULE__), "data", "model", "metabolic_reactions.csv")),
    )

    heteros = @rsubset(df, !ismissing(:Subunit))
    @select!(heteros, :RHEA_ID, :Protein, :Stoichiometry, :Subunit, :DeltaG, :RevIndex)

    homos = @rsubset(df, ismissing(:Subunit))
    @select!(homos, :RHEA_ID, :Protein, :Stoichiometry, :DeltaG, :RevIndex)

    ghomos = groupby(homos, [:RHEA_ID, :Protein])
    gheteros = groupby(heteros, [:RHEA_ID, :Subunit])

    # Build model

    model = Model()

    extend_model!(model, ghomos)
    extend_model!(model, gheteros)
    gapfill!(model)
    curate!(model) # needs to happen here, as metabolites get added that are necessary later

    add_exchanges!(model)
    add_periplasm_transporters!(model)
    add_membrane_transporters!(model)
    add_electron_transport_chain!(model)
    add_salt_transducers!(model)

    add_atpm!(model)
    add_biomass!(model)

    set_default_exchanges!(model)
    name_reactions!(model)

    # rename to conform with SBML model rules
    rename_gene_ids!(model)
    rename_metabolite_ids!(model)
    rename_reaction_ids!(model)

    model
end
