
function printmodel(model, rxns = [])

    df = DataFrame(
        rid = String[],
        Name = String[],
        Stoichiometry = String[],
        dG = Float64[],
        EC = String[],
    )

    isempty(rxns) && (rxns = collect(keys(model.reactions)))

    _f(rid) = begin
        x = model.reactions[rid].dg
        isnothing(x) ? missing : x
    end
    _g(rid) = begin
        haskey(model.reactions[rid].annotations, "KEGG_REACTION") || return ""

        x = first(model.reactions[rid].annotations["KEGG_REACTION"])
        lb = model.reactions[rid].lower_bound
        ub = model.reactions[rid].upper_bound
        
        if lb == 0 && ub != 0
            replace(x, "<=>" => "->")
        elseif lb != 0 && ub == 0
            replace(x, "<=>" => "<-")
        else
            replace(x, "<=>" => "<->")
        end
    end
    _h(rid) =  begin
        haskey(model.reactions[rid].annotations, "EC") || return ""
        ecs = model.reactions[rid].annotations["EC"]
        isnothing(ecs) && return ""
        join(ecs, "/")
    end
    _k(rid) = begin
        if startswith(rid, "EX_")
            cid = last(split(rid, "_"))
            nm = model.metabolites[cid].name
        else
            nm = model.reactions[rid].name
        end
        isnothing(nm) ? "" : nm    
    end

    for rid in rxns
        rid in keys(model.reactions) || continue 

        push!(
            df, 
            (
                rid, 
                _k(rid),
                _g(rid),
                _f(rid),
                _h(rid),
            ); 
            promote=true,
        )
    end
    
    df

    CSV.write("model.csv", df)
end
