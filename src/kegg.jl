
const EMPTYHEADER = "            " # 12 spaces for headers

get_header(h) = first(split(h, " "))

function parse_flatfile(lns)
    d = Dict{String, Vector{String}}()
    header = ""
    for ln in lns
        startswith(ln, "///") && break
        h = ln[1:12] # temp header
        header = h == EMPTYHEADER ? header : get_header(h)
        push!(get!(d, header, String[]), ln[13:end])
    end
    return d
end

function parse_reaction_stoichiometry(ln)
    stoich = Dict{String,Int64}()
    subs_prods = split(strip(ln), "<=>")
    subs = [strip(x) for x in split(subs_prods[1], " + ")]
    prods = [strip(x) for x in split(subs_prods[2], " + ")]
    for s in subs
        if startswith(s, "C") || startswith(s, "G") 
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

"""
$(TYPEDSIGNATURES)

Get the reaction name, stoichiometry, and database cross references
"""
function get_kegg_reaction(rxn_id::String; cache=true, force=false)

    !force && _is_cached("reactions", rxn_id) && return _get_cache("reactions",rxn_id)
     
    req = HTTP.request("GET", "https://rest.kegg.jp/get/$rxn_id")    
    req.status != 200 && return nothing 

    entry = parse_flatfile(split(String(req.body), "\n"))
                
    stoich = parse_reaction_stoichiometry(first(entry["EQUATION"]))

    rxn = KEGGReaction(;
            id = rxn_id,
            name = first(get(entry, "NAME", [nothing])),
            stoichiometry = stoich,
            string_stoichiometry = first(get(entry, "DEFINITION", [nothing])),
            ec = get(entry, "ENZYME", nothing),
            pathway = nothing,
            dblinks = nothing,
        )
    
    cache && _cache("reactions", rxn_id, rxn)
    
    return rxn
end

"""
$(TYPEDSIGNATURES)

Get the name and formula of a compound `met_id`.
"""
function get_kegg_compound(met_id::String; cache=true, force=false)

    !force && _is_cached("compounds", met_id) && return _get_cache("compounds", met_id)
    
    req = HTTP.request("GET", "https://rest.kegg.jp/get/$met_id")    
    req.status != 200 && return nothing 

    entry = parse_flatfile(split(String(req.body), "\n"))
    nm = first(get(entry, "NAME", [nothing]))
    nm = isnothing(nm) ? nm : replace(nm,";"=>"")  
    met = KEGGCompound(
        id = met_id,
        name = nm,
        formula = first(get(entry, "FORMULA", [nothing])),
        dblinks = nothing,
    )

    cache && _cache("compounds", met_id, met)

    return met
end

"""
$(TYPEDSIGNATURES)

Get the orthology of `ko_id`.
"""
function get_kegg_orthology(ko_id::String; cache=true, force=false)

    !force && _is_cached("orthologies", ko_id) && return _get_cache("orthologies", ko_id)
    
    req = HTTP.request("GET", "https://rest.kegg.jp/get/$ko_id")            
    req.status != 200 && return nothing

    entry = parse_flatfile(split(String(req.body), "\n"))

    ko = KEGGOrthology(
        id = ko_id,
        symbol = first(get(entry, "SYMBOL", [nothing])),
        name = replace(first(get(entry, "NAME", [nothing])),";"=>""),
        reactions = get(entry, "REACTION", nothing)
    )

    cache && _cache("orthologies", ko_id, ko)

    return ko
end
