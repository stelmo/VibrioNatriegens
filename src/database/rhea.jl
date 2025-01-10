
@kwdef mutable struct RheaReaction
    entry::String
    definition::String
    equation::Dict{String, Float64}
    enzyme::String
end

function parse_rhea_reactions(floc)
    rheas = Dict{String, RheaReaction}()
    open(floc) do IO
        entry = ""
        equation = ""
        enzyme = ""
        definition = ""
        for ln in eachline(IO)
            occursin("ENTRY", ln) && begin
                entry = string(last(split(ln, "ENTRY       RHEA:")))
            end
            occursin("DEFINITION", ln) && begin
                definition = string(last(split(ln, "DEFINITION  "))) 
            end
            occursin("EQUATION", ln) && begin
                equation = parse_reaction_stoichiometry(string(last(split(ln, "EQUATION    "))))
            end
            occursin("ENZYME", ln) && begin
                enzyme = string(last(split(ln, "ENZYME      ")))    
            end
            occursin("///", ln) && begin
                rheas[entry] = RheaReaction(; entry, equation, enzyme, definition)
                entry, equation, enzyme, definition = "", Dict{String, Float64}(), "", ""    
            end
        end
    end
    _cache("reactions", "RHEA", rheas)
end

