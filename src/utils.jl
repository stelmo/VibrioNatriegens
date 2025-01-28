
function parse_formula(x::Union{Nothing, String})
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
            name = isnothing(rxn.name) ? nothing : rxn.name

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
        haskey(model.genes, g) || begin
            model.genes[g] = Gene(name = g)            
        end
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


