
"""
$(TYPEDEF)

A struct for storing KEGG reaction information. Does not store the metabolite
information. 

$(FIELDS)
"""
@kwdef mutable struct KEGGReaction
    id::String
    name::Union{String,Nothing}
    stoichiometry::Union{Dict{String,Int64},Nothing}
    string_stoichiometry::Union{String,Nothing}
    ec::Union{Vector{String},Nothing} # multiple ECs can be assigned to a single reaction
    pathway::Union{Vector{String},Nothing} # multiple pathways possible
    dblinks::Union{Vector{String},Nothing}
end

"""
$(TYPEDEF)

A struct for storing KEGG compound information.

$(FIELDS)
"""
@kwdef mutable struct KEGGCompound
    id::String
    name::Union{String,Nothing}
    molarmass::Union{Float64,Nothing}
    formula::Union{String,Nothing}
end

"""
$(TYPEDEF)

A struct for storing KEGG orthology information.

$(FIELDS)
"""
@kwdef mutable struct KEGGOrthology
    id::String
    symbol::String
    name::String
    reactions::Union{Vector{String},Nothing}
end

const EMPTYHEADER = "            " # 12 spaces for headers

get_header(h) = first(split(h, " "))

function parse_flatfile(lns)
    d = Dict{String,Vector{String}}()
    header = ""
    for ln in lns
        startswith(ln, "///") && break
        startswith(ln, "REFERENCE") && continue
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

Get kegg entry.
"""
function get_kegg_entry(id::String, db::String; cache = true, force = false)
    if !force && _is_cached(db, id)
        entry = _get_cache(db, id)
    else
        sleep_ms(333) # from Timers.jl, KEGG allows max of 3 queries per second
        req = HTTP.request("GET", "https://rest.kegg.jp/get/$id")
        req.status != 200 && return nothing
        entry = parse_flatfile(split(String(req.body), "\n"))
        cache && _cache(db, id, entry)
    end
    entry
end

"""
$(TYPEDSIGNATURES)

Get the reaction name, stoichiometry, and database cross references
"""
function get_kegg_reaction(rxn_id::String; cache = true, force = false)

    entry = VibrioNatriegens.get_kegg_entry(rxn_id, "reactions"; cache, force)
    name = first(get(entry, "NAME", [nothing]))
    name = isnothing(name) ? nothing : replace(name, ";" => "")
    stoichiometry = VibrioNatriegens.parse_reaction_stoichiometry(first(entry["EQUATION"]))
    pathway = haskey(entry, "PATHWAY") ? [x[3:7] for x in entry["PATHWAY"]] : nothing
    dblinks = haskey(entry, "DBLINKS") ? [x for x in entry["DBLINKS"]] : nothing
    ec = get(entry, "ENZYME", nothing)
    ec = isnothing(ec) ? nothing : filter(!isempty, split(first(ec), " "))
    string_stoichiometry = first(get(entry, "DEFINITION", [nothing]))

    rxn = KEGGReaction(;
        id = rxn_id,
        name,
        stoichiometry,
        string_stoichiometry,
        ec,
        pathway,
        dblinks,
    )

    return rxn
end

"""
$(TYPEDSIGNATURES)

Get the name and formula of a compound `met_id`.
"""
function get_kegg_compound(met_id::String; cache = true, force = false)

    entry = VibrioNatriegens.get_kegg_entry(met_id, "compounds"; cache, force)

    name = first(get(entry, "NAME", [nothing]))
    name = isnothing(name) ? name : replace(name, ";" => "")
    molarmass = first(get(entry, "MOL_WEIGHT", [nothing]))
    molarmass = isnothing(molarmass) ? nothing : parse(Float64, molarmass)

    met = KEGGCompound(;
        id = met_id,
        name,
        molarmass,
        formula = first(get(entry, "FORMULA", [nothing])),
    )

    return met
end

"""
$(TYPEDSIGNATURES)

Get the orthology of `ko_id`.
"""
function get_kegg_orthology(ko_id::String; cache = true, force = false)

    entry = get_kegg_entry(ko_id, "orthologies"; cache, force)

    ko = KEGGOrthology(
        id = ko_id,
        symbol = first(get(entry, "SYMBOL", [nothing])),
        name = replace(first(get(entry, "NAME", [nothing])), ";" => ""),
        reactions = get(entry, "REACTION", nothing),
    )

    return ko
end

"""
$(TYPEDSIGNATURES)

List all modules in reaction.
"""
function list_kegg_modules()
    req = HTTP.request("GET", "https://rest.kegg.jp/list/module")
    req.status != 200 && return nothing
    Dict(
        string(first(split(x, "\t"))) => string(last(split(x, "\t"))) for
        x in split(String(req.body), "\n") if x != ""
    )
end

"""
$(TYPEDSIGNATURES)

Get module info.
"""
get_kegg_module(modid; cache = true, force = false) =
    get_kegg_entry(modid, "modules"; cache, force)

"""
$(TYPEDSIGNATURES)

"""
function get_kegg_reactions_in_pathway(mapid; cache = true, force = false)
    if !force && _is_cached("link_path_rn", mapid)
        entry = _get_cache("link_path_rn", mapid)
    else
        req = HTTP.request("GET", "https://rest.kegg.jp/link/rn/$mapid")
        req.status != 200 && return nothing
        entry = [string(x[end-5:end]) for x in split(String(req.body), "\n") if x != ""]
        cache && _cache("link_path_rn", mapid, entry)
    end
    entry
end

"""
$(TYPEDSIGNATURES)
"""
function parse_kegg_definition(definition_string, seperator = ' ')
    terms = String[]
    bracket_level = 0
    current_term = ""
    for i = 1:length(definition_string)

        bracket_level += definition_string[i] == '(' ? 1 : 0
        bracket_level -= definition_string[i] == ')' ? 1 : 0

        if bracket_level != 0 || (bracket_level == 0 && definition_string[i] != seperator)
            current_term *= definition_string[i]
        end

        if bracket_level == 0 &&
           (definition_string[i] == seperator || i == length(definition_string))
            push!(terms, current_term)
            current_term = ""
        end
    end
    terms
end

function get_kegg_reactions_from_map(mapid; cache = true, force = false)
    mapid = String(mapid)
    if !force && _is_cached("maps", mapid)
        entry = _get_cache("maps", mapid)
    else
        sleep_ms(333) # from Timers.jl, KEGG allows max of 3 queries per second
        req = HTTP.request("GET", "https://rest.kegg.jp/link/rn/$mapid")
        req.status != 200 && return nothing
        entry = String[]
        for ln in split(String(req.body), "\n")
            en = last(split(ln, "\t"))
            en != "" && push!(entry, en[4:end])
        end
        entry = unique(entry)
        cache && _cache("maps", mapid, entry)
    end
    entry
end
