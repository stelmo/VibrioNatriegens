
"""
$(TYPEDSIGNATURES)

Build the model.
"""
function build_model()

    df = DataFrame(CSV.File(joinpath("data", "model", "metabolic_reactions.csv")))

    heteros = @rsubset(df, !ismissing(:Subunit))
    @select!(heteros, :RHEA_ID, :Protein, :Stoichiometry, :Subunit, :DeltaG, :RevIndex)
    
    homos = @rsubset(df, ismissing(:Subunit))
    @select!(homos, :RHEA_ID, :Protein, :Stoichiometry, :DeltaG, :RevIndex)

    ghomos = groupby(homos, [:RHEA_ID, :Protein])
    gheteros = groupby(heteros, [:RHEA_ID, :Subunit])

    # Build model

    model = Model()

    VibrioNatriegens.extend_model!(model, ghomos)
    VibrioNatriegens.extend_model!(model, gheteros)
    VibrioNatriegens.gapfill!(model)
    VibrioNatriegens.curate!(model) # needs to happen here, as metabolites get added that are necessary later

    VibrioNatriegens.add_exchanges!(model)
    VibrioNatriegens.add_periplasm_transporters!(model)
    VibrioNatriegens.add_membrane_transporters!(model)
    VibrioNatriegens.add_electron_transport_chain!(model)
    VibrioNatriegens.add_salt_transducers!(model)

    VibrioNatriegens.add_atpm!(model)
    VibrioNatriegens.add_biomass!(model)

    VibrioNatriegens.set_default_exchanges!(model)
    VibrioNatriegens.name_reactions!(model)
    
    model
end
