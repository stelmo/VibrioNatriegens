module VibrioNatriegens

using HTTP, DocStringExtensions, Scratch, Serialization

const CACHE_DIRS = ["reactions", "compounds", "orthologies"]
CACHE_LOCATION = ""

function __init__()
    global CACHE_LOCATION = @get_scratch!("kegg_data")

    for dir in CACHE_DIRS
        !isdir(joinpath(CACHE_LOCATION, dir)) && mkdir(joinpath(CACHE_LOCATION, dir))
    end

end

include("cache.jl")
include("types.jl")
include("kegg.jl")

using CSV, AbstractFBCModels, JSONFBCModels, JSON
import AbstractFBCModels as A

include("auto.jl")
include("exchanges.jl")
include("build_model.jl")

export build_model

end
