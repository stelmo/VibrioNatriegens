using JSON
using VibrioNatriegens
using AbstractFBCModels
import AbstractFBCModels as A

model = VibrioNatriegens.build_model()
prefixes = ["EX_", "BIOMASS", "ATPM", "DF_", "SYM_", "ANTI_", "PERM_", "ABC_", "PTS"]
metabolic_rxns = filter(x -> !any(startswith.(x, prefixes)), A.reactions(model))


map_files = [
    "aminoacids.json",
    "carbohydrates.json",
    "cofactors.json",
    "energy.json",
    "lipids map.json",
    "lipids.json",
    "nucleotides.json",
]
map_bigg_ids = Set{String}()
for file in map_files
    println("Reading $file")
    eschermap = JSON.parsefile(joinpath("maps", file))[2]
    for (_, rxn) in eschermap["reactions"]
        push!(map_bigg_ids, rxn["bigg_id"])
    end
end

unmapped_rxns = setdiff(metabolic_rxns, map_bigg_ids)

open("unmapped_rxns.txt", "w") do io
    write(io, join(unmapped_rxns, "\n"))
end

println("Results written to rhea_unmapped.txt")
