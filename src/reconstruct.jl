
"""
$(TYPEDSIGNATURES)

Build the model.
"""
function build_model()

    df = DataFrame(
        CSV.File(joinpath(pkgdir(@__MODULE__), "data", "model", "metabolic_reactions.csv")),
    )

    heteros = @rsubset(df, !ismissing(:Isozyme))
    @select!(heteros, :RHEA_ID, :Protein, :Stoichiometry, :Isozyme, :DeltaG, :RevIndex)

    homos = @rsubset(df, ismissing(:Isozyme))
    @select!(homos, :RHEA_ID, :Protein, :Stoichiometry, :DeltaG, :RevIndex)

    ghomos = groupby(homos, [:RHEA_ID, :Protein])
    gheteros = groupby(heteros, [:RHEA_ID, :Isozyme])

    # Build model

    model = Model()

    extend_model!(model, ghomos)
    extend_model!(model, gheteros)
    gapfill!(model)
    curate!(model) # needs to happen here, as metabolites get added that are necessary later

    # these are all manually added reactions
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

    fix_noncytosolic_metabolite_annotations!(model)
    set_default_exchanges!(model)
    name_reactions_genes!(model)

    # rename to conform with SBML model rules
    rename_gene_ids!(model)
    rename_metabolite_ids!(model)
    rename_reaction_ids!(model)

    model
end
