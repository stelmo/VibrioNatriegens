
"""
$(TYPEDSIGNATURES)

Parse a string metabolite formula into a dictionary.
"""
function parse_formula(x::A.Maybe{String})
    isnothing(x) && return nothing
    x == "" && return nothing

    res = Dict{String,Int}()
    pattern = @r_str "([A-Z][a-z]*)([1-9][0-9]*)?"
    for m in eachmatch(pattern, x)
        res[m.captures[1]] = isnothing(m.captures[2]) ? 1 : parse(Int, m.captures[2])
    end
    return res
end

"""
$(TYPEDSIGNATURES)

Add metabolic reactions to the model.
"""
function extend_model!(model, row)

    rid = string(row.Rhea) # rhea ID, required
    subunit_stoichiometry = row.Stoichiometry # gene stoichiometry in isozyme, can be missing
    protein = row.Protein # gene ID, can be missing
    isozyme = row.Isozyme # the isozyme ID, can be missing

    if haskey(model.reactions, rid) # extend a current reaction (only isozymes can get added at this point)
        @assert !ismissing(protein) "Protein ID is missing for reaction $rid"
        @assert !ismissing(subunit_stoichiometry) "Subunit stoichiometry is missing for reaction $rid"
        @assert !ismissing(isozyme) "Isozyme ID is missing for reaction $rid"

        isos = model.reactions[rid].gene_association
        if haskey(isos, isozyme) # add to complex
            isos[isozyme].gene_product_stoichiometry[protein] = subunit_stoichiometry
        else # new isozyme
            isos[isozyme] =
                Isozyme(Dict(protein => subunit_stoichiometry), nothing, nothing)
        end

        add_gene!(model, protein) # add gene to model
    else # add a new reaction
        rxn = get_reaction(rid) # rhea lookup, from cache
        met_stoichiometry = rxn.stoichiometry

        iso =
            ismissing(protein) ? nothing :
            Dict(
                isozyme =>
                    Isozyme(Dict(protein => subunit_stoichiometry), nothing, nothing),
            )

        model.reactions[rid] = Reaction(;
            lower_bound = -1000.0, # reversible by default
            upper_bound = 1000.0,
            gene_association = iso,
            stoichiometry = met_stoichiometry,
            annotations = Dict("rhea-reaction-description" => [rxn.equation]),
        )
        add_metabolites!(model, collect(keys(met_stoichiometry))) # add metabolites to model
        ismissing(protein) || add_gene!(model, protein) # add gene to model if it is available, other gap
    end

end

add_gene!(model, gid) = haskey(model.genes, gid) || (model.genes[gid] = Gene(;))

add_genes!(model, gids) = begin
    for gid in gids
        add_gene!(model, gid)
    end
end

function add_metabolites!(model, ms)
    for m in ms
        haskey(model.metabolites, m) || begin
            met = get_metabolite(m)
            model.metabolites[m] = Metabolite(;
                name = met.name,
                formula = parse_formula(met.formula),
                charge = met.charge,
                compartment = "Cytosol",
            )
        end
    end
end
