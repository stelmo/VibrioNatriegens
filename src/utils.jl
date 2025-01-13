
function parse_formula(x::Union{Nothing, String})
    isnothing(x) && return nothing
    x == "" && return nothing

    res = Dict{String,Int}()
    pattern = @r_str "([A-Z][a-z]*)([1-9][0-9]*)?"
    for m in eachmatch(pattern, x)
        res[m.captures[1]] = isnothing(m.captures[2]) ? 1 : parse(Int, m.captures[2])
    end
    return res
end
