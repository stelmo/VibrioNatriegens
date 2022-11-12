ecs = String[]
open(joinpath("data", "deepec", "DeepEC_Result.txt")) do io
    firstline =  true
    for ln in eachline(io)
        firstline && (firstline = false; continue)
        prts = split(string(last(split(ln, "\t"))), ":")
        push!(ecs, last(prts))
    end
end
unique(ecs)

open(joinpath("data", "deepec", "only_ecs.txt"), "w") do io
    for ec in ecs
        write(io, ec, "\n")
    end
end
