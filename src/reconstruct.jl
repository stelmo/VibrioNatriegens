
"""
$(TYPEDSIGNATURES)

Build the genome-scale metabolic model of Vibrio natriegens from scratch. This
is usually time intensive, so it is recommended to build the model once, and
then serialize it for quick access.
"""
function build_model()
    #=
    It is very important that every reaction listed in metabolic_reactions.csv
    is the reference RHEA reaction. This typically corresponds to the smallest
    number in the RHEA quartet.
    
    =#

    # add reactions with genetic evidence
    core_rows = CSV.File(joinpath(pkgdir(@__MODULE__), "data", "model", "reactions_metabolic.csv"))

    # add gap filled reactions    
    gap_rows = CSV.File(joinpath(pkgdir(@__MODULE__), "data", "model", "reactions_gaps.csv"))

    # per-cache everything
    get_reactions(core_rows.Rhea)
    get_reactions(gap_rows.Rhea)

    # Build model
    # model = Model()

    model = Model()

    for row in [core_rows; gap_rows]
        extend_model!(model, row)
    end

    # curate the model further
    curate!(model) # needs to happen here, as metabolites get added that are necessary later
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

    name_genes!(model)
    name_reactions!(model)
    acronym_reactions!(model)

    set_default_exchanges!(model)
    switch_off_salt_transporters!(model)
    deactivate_antiports!(model)
    haskey(model.reactions, "PERM_29101") && delete!(model.reactions, "PERM_29101") # house keeping for Na+

    check_gene_names(model)
    check_reaction_acronyms(model)
    check_rhea_ref(model)

    model
end
