function guess_rxn_split(ln)
    occursin("<=>", ln) && return "<=>"
    occursin("=>", ln) && return "=>"
    occursin("<=", ln) && return "<="
    occursin("=", ln) && return "="
    nothing
end

"""
$(TYPEDSIGNATURES)

Parse a reaction stoichiometry string.
"""
function parse_reaction_stoichiometry(ln)
    rxn_eq = guess_rxn_split(ln)

    stoich = Dict{String,Int64}()
    subs_prods = split(strip(ln), rxn_eq)
    subs = [strip(x) for x in split(subs_prods[1], " + ")]
    prods = [strip(x) for x in split(subs_prods[2], " + ")]
    for s in subs
        if startswith(s, "C") || startswith(s, "G") # happily also works for chebi/rhea
            id = s
            v = 1
        else
            x = split(s)
            id = String(x[2])
            v = isnothing(tryparse(Int64, x[1])) ? 0 : tryparse(Int64, x[1])
        end
        stoich[id] = get!(stoich, id, 0) - v
    end

    for p in prods
        if startswith(p, "C") || startswith(p, "G")
            id = p
            v = 1
        else
            x = split(p)
            id = String(x[2])
            v = isnothing(tryparse(Int64, x[1])) ? 0 : tryparse(Int64, x[1])
        end
        stoich[id] = get!(stoich, id, 0) + v
    end
    stoich
end
