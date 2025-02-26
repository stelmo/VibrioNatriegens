
function print_reactions(model::Model, rids::Vector{String}, file_name)
    df = DataFrame(
        rid = String[],
        Name = String[],
        Stoichiometry = String[],
        dG = Float64[],
        EC = String[],
    )

    _dg(rid) = begin
        x = model.reactions[rid].dg
        isnothing(x) ? missing : x
    end

    _stoichiometry(rid) = begin

        ss = model.reactions[rid].stoichiometry

        s1 = join(
            [string(v) * "*" * A.metabolite_name(model, k) for (k, v) in ss if v < 0],
            " + ",
        )
        s2 = join(
            [string(v) * "*" * A.metabolite_name(model, k) for (k, v) in ss if v > 0],
            " + ",
        )

        lb = model.reactions[rid].lower_bound
        ub = model.reactions[rid].upper_bound

        if lb == 0 && ub != 0
            dir = " -> "
        elseif lb != 0 && ub == 0
            dir = " <- "
        else
            dir = " <-> "
        end

        s1 * dir * s2
    end

    _ec(rid) = begin
        haskey(model.reactions[rid].annotations, "EC") || return ""
        ecs = model.reactions[rid].annotations["EC"]
        isnothing(ecs) && return ""
        join(last.(split.(ecs, "/")), "/")
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

    for rid in rids
        # println(rid)
        push!(
            df,
            (rid, _name(rid), _stoichiometry(rid), _dg(rid), _ec(rid));
            promote = true,
        )
    end

    df

    CSV.write(file_name, df)
end

print_reactions(model::Model) =
    print_reactions(model, A.reactions(model), "reactions-model.csv")

function print_metabolites(model)

    mids = A.metabolites(model)
    charges = Union{Missing,Int}[]
    formulas = Union{Missing,String}[]
    names = Union{Missing,String}[]

    getx(x) = isnothing(x) ? missing : x

    for mid in mids
        push!(charges, getx(A.metabolite_charge(model, mid)))
        push!(names, getx(A.metabolite_name(model, mid)))

        f = A.metabolite_formula(model, mid)

        push!(formulas, isnothing(f) ? missing : join(k * string(v) for (k, v) in f))
    end

    df = DataFrame(Metabolite = mids, Name = names, Charge = charges, Formula = formulas)
    CSV.write("metabolites-model.csv", df)
end

