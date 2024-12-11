using DataFrames, CSV, DataFramesMeta

p = joinpath("data", "annotations", "hamap", "trusted.txt")

s = read(p, String)

entries = []
gid = ""
ss = ""
subunitstarted = false
for ln in eachline(p)
    if startswith(ln, "//")
        if ss != ""
            push!(entries, (; gene = gid, subunit = ss))
        end
        gid = ""
        ss = ""
    end
    if startswith(ln, "CC   -!- SUBUNIT:")
        ss *= last(split(ln, "-!- SUBUNIT: "))
        subunitstarted = true
    elseif subunitstarted && !startswith(ln, "CC   -!-")
        println(ln)
        ss *= last(split(ln, "C       "))
    else
        subunitstarted = false
    end
    if startswith(ln, "**HA Submitted Name: ")
        gid = last(split(ln, "Submitted Name: "))
    end
end


df = DataFrame(entries)
@rename!(df, :Protein = :gene, :Subunit = :subunit)
CSV.write(joinpath("data", "annotations", "hamap", "hamap_subunits.csv"), df)
