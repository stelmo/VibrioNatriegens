"""
$(TYPEDSIGNATURES)

Add metabolic reactions to the model.
"""
function extend_model!(model, dfs)

    gs = String[]
    ms = String[]

    for df in dfs
        rid = first(df.Reaction)
        grr = String.(df.Protein[:])
        stoich = Int.(df.Stoichiometry[:])
        append!(gs, grr)


        iso = X.Isozyme(; gene_product_stoichiometry = Dict(grr .=> stoich))

        if haskey(model.reactions, rid) # isozyme
            push!(model.reactions[rid].gene_association, iso)

        else # first time seeing this reaction
            kegg_rxn = get_kegg_reaction(rid)

            _ms = string.(keys(kegg_rxn.stoichiometry))
            for m in _ms
                occursin("(", m) && println(rid)
            end

            append!(ms, _ms)
            ecs = isnothing(kegg_rxn.ec) ? [""] : kegg_rxn.ec

            model.reactions[rid] = Reaction(;
                name = kegg_rxn.name,
                lower_bound = -1000,
                upper_bound = 1000,
                gene_association = [iso],
                stoichiometry = kegg_rxn.stoichiometry,
                annotations = Dict(
                    "KEGG_REACTION" => [kegg_rxn.string_stoichiometry],
                    "EC" => ecs,
                ),
            )
        end
    end

    # add genes
    for g in unique(gs)
        model.genes[g] = Gene(name = g)
    end

    # add metabolites
    for m in unique(ms)
        kegg_met = get_kegg_compound(m)
        model.metabolites[m] = Metabolite(
            name = kegg_met.name,
            formula = parse_formula(kegg_met.formula),
            compartment = "Cytosol",
        )
    end

end

"""
$(TYPEDSIGNATURES)

Build the model.
"""
function build_model()

    df = DataFrame(
        XLSX.readtable(
            joinpath("data", "curation", "curated", "base_reactions.xlsx"),
            "metabolism",
        ),
    )

    heteros = @rsubset(df, !ismissing(:Subunit))
    homos = @rsubset(df, ismissing(:Subunit))
    @select!(homos, :Reaction, :Protein, :Stoichiometry)

    ghomos = groupby(homos, [:Reaction, :Protein])
    gheteros = groupby(heteros, [:Reaction, :Subunit])

    #! Build model

    model = VibrioNatriegens.Model()

    VibrioNatriegens.extend_model!(model, ghomos)
    VibrioNatriegens.extend_model!(model, gheteros)
    VibrioNatriegens.reaction_directions!(model) # only metabolic reactions
    VibrioNatriegens.curate!(model)

    VibrioNatriegens.add_exchanges!(model)
    VibrioNatriegens.add_periplasm_transporters!(model)
    VibrioNatriegens.add_membrane_transporters!(model)
    VibrioNatriegens.add_electron_transport_chain!(model)
    VibrioNatriegens.add_salt_transducers!(model)

    VibrioNatriegens.add_atpm!(model)
    VibrioNatriegens.add_biomass!(model)

    VibrioNatriegens.printmodel(model)
    model
end
