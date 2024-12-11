module VibrioNatriegens

using HTTP, DocStringExtensions, Scratch, Serialization, Timers
import AbstractFBCModels as A
import JSONFBCModels: parse_formula
import SparseArrays: sparse, findnz
import COBREXA as X
using eQuilibrator, Measurements, Unitful
using DataFrames, CSV


const CACHE_DIRS =
    ["reactions", "compounds", "orthologies", "modules", "maps", "directionality","link_path_rn"]
CACHE_LOCATION = ""

function __init__()
    global CACHE_LOCATION = @get_scratch!("kegg_data")
    for dir in CACHE_DIRS
        !isdir(joinpath(CACHE_LOCATION, dir)) && mkdir(joinpath(CACHE_LOCATION, dir))
    end
end

include("cache.jl")
include("kegg.jl")

include("model.jl")
include("build.jl")
include("curate.jl")
include("biomass.jl")
include("exchanges.jl")
include("directions.jl")
include("printmodel.jl")

end
