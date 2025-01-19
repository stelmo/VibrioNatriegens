
"""
$(TYPEDSIGNATURES)

Build the model.
"""
function build_model()

    df = DataFrame(CSV.File(joinpath("data", "model", "metabolic_reactions.csv")))

    heteros = @rsubset(df, !ismissing(:Subunit))
    @select!(heteros, :RHEA_ID, :Protein, :Stoichiometry, :Subunit)
    
    homos = @rsubset(df, ismissing(:Subunit))
    @select!(homos, :RHEA_ID, :Protein, :Stoichiometry,)

    ghomos = groupby(homos, [:RHEA_ID, :Protein])
    gheteros = groupby(heteros, [:RHEA_ID, :Subunit])

    # Build model

    model = Model()

    VibrioNatriegens.extend_model!(model, ghomos)
    VibrioNatriegens.extend_model!(model, gheteros)
    # VibrioNatriegens.reaction_directions!(model) # only metabolic reactions
    # VibrioNatriegens.curate!(model)

    # VibrioNatriegens.add_exchanges!(model)
    # VibrioNatriegens.add_periplasm_transporters!(model)
    # VibrioNatriegens.add_membrane_transporters!(model)
    # VibrioNatriegens.add_electron_transport_chain!(model)
    # VibrioNatriegens.add_salt_transducers!(model)

    # VibrioNatriegens.add_atpm!(model)
    # VibrioNatriegens.add_biomass!(model)

    # VibrioNatriegens.printmodel(model)
    model
end

"""
$(TYPEDSIGNATURES)

Add metabolic reactions to the model.
"""
function extend_model!(model, dfs)

    gs = String[]
    ms = RheaReactions.RheaMetabolite[]

    for df in dfs

        rid = first(df.RHEA_ID)
        grr = String.(df.Protein[:])
        stoich = Int.(df.Stoichiometry[:])
        append!(gs, grr)

        iso = X.Isozyme(; gene_product_stoichiometry = Dict(grr .=> stoich))

        if haskey(model.reactions, string(rid)) # isozyme

            push!(model.reactions[string(rid)].gene_association, iso)

        else # first time seeing this reaction
            
            rxn = get_reaction(rid)

            coeff_mets = get_reaction_metabolites(rid)
            stoichiometry = Dict(
                string(v.accession) => s
                for (s, v) in coeff_mets
            )

            append!(ms, last.(coeff_mets))
            
            ecs = isnothing(rxn.ec) ? [""] : rxn.ec
            name = isnothing(rxn.name) ? "" : rxn.name

            model.reactions[string(rid)] = Reaction(;
                name = name,
                lower_bound = -1000,
                upper_bound = 1000,
                gene_association = [iso],
                stoichiometry = stoichiometry,
                annotations = Dict(
                    "REACTION" => [rxn.equation],
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
    for m in ms
        haskey(model.metabolites, m.accession) || begin
                model.metabolites[m.accession] = Metabolite(;
                name = m.name,
                formula = parse_formula(m.formula),
                charge = m.charge,
                compartment = "Cytosol",
            )
        end
    end
end


