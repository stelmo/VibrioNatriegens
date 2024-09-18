"""
$(TYPEDSIGNATURES)

Build the genome-scale metabolic model.
"""
function build_model()

    # create draft based on KEGG's annotations
    model = VibrioNatriegens.automatic_build()

    # add exchanges
    VibrioNatriegens.add_exchanges!(model)

end
