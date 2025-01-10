module VibrioNatriegens

using HTTP, DocStringExtensions, Scratch, Serialization, Timers
import AbstractFBCModels as A
import JSONFBCModels: parse_formula
import SparseArrays: sparse, findnz
import COBREXA as X
using eQuilibrator, Measurements, Unitful
using DataFrames, CSV, DataFramesMeta, XLSX
import ConstraintTrees as C

const CACHE_DIRS = [
    "reactions",
    "compounds",
    "orthologies",
    "modules",
    "maps",
    "directionality",
    "link_path_rn",
]
CACHE_LOCATION = ""

function __init__()
    global CACHE_LOCATION = @get_scratch!("kegg_data")
    for dir in CACHE_DIRS
        !isdir(joinpath(CACHE_LOCATION, dir)) && mkdir(joinpath(CACHE_LOCATION, dir))
    end
end

include("database/cache.jl")
include("database/kegg.jl")
include("database/utils.jl")
include("database/rhea.jl")
include("database/chebi.jl")


include("model/model.jl")
include("model/printmodel.jl")
include("model/build.jl")
include("model/biomass.jl")
include("model/transporters.jl")
include("model/directions.jl")
include("model/curate.jl")

end
