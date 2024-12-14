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

            append!(ms, string.(keys(kegg_rxn.stoichiometry)))
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
        if occursin("(", m)
            println(m)
        else
            kegg_met = get_kegg_compound(m)
            model.metabolites[m] =
                Metabolite(name = kegg_met.name, formula = parse_formula(kegg_met.formula))
        end
    end

end

function build_model(; pathway_maps = [], add_rxns = [])

    df = DataFrame(
        XLSX.readtable(
            joinpath("data", "curation", "curated", "base_reactions.xlsx"),
            "Sheet1",
        ),
    )

    # this is useful for troubleshooting
    rxns = String[]
    for m in pathway_maps
        append!(rxns, get_kegg_reactions_in_pathway(m))
    end
    append!(rxns, add_rxns)

    !isempty(rxns) && @rsubset!(df, :Reaction in rxns)

    heteros = @rsubset(df, !ismissing(:Subunit))
    homos = @rsubset(df, ismissing(:Subunit))
    @select!(homos, :Reaction, :Protein, :Stoichiometry)

    ghomos = groupby(homos, [:Reaction, :Protein])
    gheteros = groupby(heteros, [:Reaction, :Subunit])

    #! Build model

    model = Model()

    extend_model!(model, ghomos)
    extend_model!(model, gheteros)
    reaction_directions!(model)
    add_exchanges!(model)
    add_atpm!(model)
    add_biomass!(model)
    curate!(model)

    printmodel(model)
    model
end
