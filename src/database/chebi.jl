
@kwdef mutable struct CHEBICompound
    entry::String
    mass::X.Maybe{Float64} = nothing
    charge::X.Maybe{Int64} = nothing
    formula::X.Maybe{String} = nothing
    name::X.Maybe{String} = nothing
    reactions::X.Maybe{String} = nothing
end

function parse_chebi_compounds(floc, floc2)
    chebis = Dict{String, CHEBICompound}()
    open(floc) do IO
        for ln in eachline(IO)
            _, entry,_,type,chem_data = string.(split(ln, "\t"))
            d = get!(chebis, entry, CHEBICompound(; entry))
            type == "FORMULA" && begin
                d.formula = chem_data
            end
            type == "MASS" && begin
                d.mass = parse(Float64, chem_data)               
            end
            type == "CHARGE" && begin
                d.charge = parse(Int, chem_data) 
            end
        end
    end

    open(floc2) do IO
        entry, name, reaction = "", "", ""
        for ln in eachline(IO)
            occursin("ENTRY", ln) && begin
                entry = string(last(split(ln, "ENTRY       CHEBI:")))
            end
            occursin("NAME", ln) && begin
                name = string(last(split(ln, "NAME        "))) 
            end
            occursin("REACTION", ln) && begin
                reaction = string(last(split(ln, "REACTION    ")))
            end
            occursin("///", ln) && begin
                if haskey(chebis, entry)
                    chebis[entry].name = name == "" ? nothing : name
                    chebis[entry].reactions = reaction == "" ? nothing : reaction                    
                end
            end
        end
    end

    _cache("compounds", "CHEBI", chebis) 
end
