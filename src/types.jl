
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
    string_stoichiometry::Union{String, Nothing}
    ec::Union{Vector{String},Nothing} # multiple ECs can be assigned to a single reaction
    pathway::Union{Vector{String},Nothing} # multiple pathways possible
    dblinks::Union{Dict{String,Vector{String}},Nothing}
end

"""
$(TYPEDEF)

A struct for storing KEGG compound information.

$(FIELDS)
"""
@kwdef mutable struct KEGGCompound
    id::String
    name::Union{String, Nothing}
    formula::Union{String,Nothing}
    dblinks::Union{Dict{String,Vector{String}},Nothing}
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
    reactions::Union{Vector{String}, Nothing}
end
