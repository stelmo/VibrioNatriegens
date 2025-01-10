
function printmodel(model::Model, rxns = [])

    df = DataFrame(
        rid = String[],
        Name = String[],
        Stoichiometry = String[],
        dG = Float64[],
        EC = String[],
    )

    isempty(rxns) && (rxns = collect(keys(model.reactions)))

    _dg(rid) = begin
        x = model.reactions[rid].dg
        isnothing(x) ? missing : x
    end
    _stoichiometry(rid) = begin

        ss = model.reactions[rid].stoichiometry

        s1 = join([A.metabolite_name(model, k) for (k, v) in ss if v < 0], " + ")
        s2 = join([A.metabolite_name(model, k) for (k, v) in ss if v > 0], " + ")
        x = s1 * " <=> " * s2

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
    _ec(rid) = begin
        haskey(model.reactions[rid].annotations, "EC") || return ""
        ecs = model.reactions[rid].annotations["EC"]
        isnothing(ecs) && return ""
        join(ecs, "/")
    end
    _name(rid) = begin
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
            (rid, _name(rid), _stoichiometry(rid), _dg(rid), _ec(rid));
            promote = true,
        )
    end

    df

    CSV.write("model.csv", df)
end
