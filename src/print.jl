
function print_reactions(model::Model, rids::Vector{String})

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
        haskey(model.reactions[rid].annotations, "rhea-ec") || return ""
        ecs = model.reactions[rid].annotations["rhea-ec"]
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

    open("vnat_rxns.tsv", "w") do io
        write(io, join(["ID", "Name", "Reaction", "EC"], "\t"), "\n")
        for rid in rids
            write(io, join([rid, _name(rid), _stoichiometry(rid), _ec(rid)], "\t"), "\n")
        end
    end
end

print_reactions(model::Model) = print_reactions(model, A.reactions(model))

function print_metabolites(model)

    getx(x) = isnothing(x) ? "" : x

    mids = A.metabolites(model)

    open("vnat_mets.csv", "w") do io
        write(io, join(["ID", "Name", "Charge", "Formula"], "\t"), "\n")
        for mid in mids
            f = A.metabolite_formula(model, mid)
            f = isnothing(f) ? "" : join(k * string(v) for (k, v) in f)
            write(
                io,
                join(
                    [
                        mid,
                        getx(A.metabolite_name(model, mid)),
                        getx(A.metabolite_charge(model, mid)),
                        f,
                    ],
                    "\t",
                ),
                "\n",
            )
        end
    end
end

