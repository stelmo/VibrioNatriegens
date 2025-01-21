
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
    VibrioNatriegens.curate!(model) # needs to happen here, as metabolites get added that are necessary later

    VibrioNatriegens.add_exchanges!(model)
    VibrioNatriegens.add_periplasm_transporters!(model)
    VibrioNatriegens.add_membrane_transporters!(model)
    VibrioNatriegens.add_electron_transport_chain!(model)
    VibrioNatriegens.add_salt_transducers!(model)

    VibrioNatriegens.add_atpm!(model)
    # VibrioNatriegens.add_biomass!(model)

    VibrioNatriegens.set_default_exchanges!(model)

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

            # direction
            reversibility_index_threshold = 3.0
            rev_ind = ismissing(first(df.RevIndex)) ? nothing : first(df.RevIndex) 
            dg = ismissing(first(df.DeltaG)) ? nothing : first(df.DeltaG) 


            if isnothing(rev_ind) || (abs(rev_ind) <= reversibility_index_threshold)
                lb = -1000
                ub = 1000
            elseif rev_ind < -reversibility_index_threshold # forward
                lb = 0
                ub = 1000
            elseif rev_ind > reversibility_index_threshold # reverse
                lb = -1000
                ub = 0
            end

            model.reactions[string(rid)] = Reaction(;
                name = name,
                lower_bound = lb,
                upper_bound = ub,
                dg = dg,
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


