
function parse_formula(x::Union{Nothing,String})
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
function extend_model!(model, dfs)

    gs = String[]
    ms = String[]

    for df in dfs

        rid = first(df.RHEA_ID)

        if any(isnothing.(df.Protein[:]))
            iso = nothing
        else
            grr = String.(df.Protein[:])
            stoich = Int.(df.Stoichiometry[:])
            append!(gs, grr)
            iso = X.Isozyme(; gene_product_stoichiometry = Dict(grr .=> stoich))
        end

        if haskey(model.reactions, string(rid)) # isozyme
            push!(model.reactions[string(rid)].gene_association, iso)

        else # first time seeing this reaction

            rxn = get_reaction(rid; verbose = false)
            stoichiometry = rxn.stoichiometry
            append!(ms, collect(keys(stoichiometry)))

            model.reactions[string(rid)] = Reaction(;
                lower_bound = -1000.0, # reversible by default, equilibrator causes more problems than it fixes
                upper_bound = 1000.0,
                gene_association = isnothing(iso) ? nothing : [iso],
                stoichiometry = stoichiometry,
                annotations = Dict("rhea-reaction-description" => [rxn.equation]),
            )
        end
    end

    add_genes!(model, gs)
    add_metabolites!(model, ms)
end

function add_genes!(model, gs)
    for g in unique(gs)
        haskey(model.genes, g) || begin
            model.genes[g] = Gene(;)
        end
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
