module VibrioNatriegens

using RheaReactions, CSV # load reactions
import AbstractFBCModels as A # model type
import SparseArrays as S # model accessors
using DocStringExtensions

include("reconstruct_utils.jl")
include("model.jl")
include("gene_annotations.jl")
include("metabolite_annotations.jl")
include("reaction_annotations.jl")
include("print.jl") # useful for troubleshooting

include("reconstruct.jl")
include("curate.jl")
include("transporters.jl")
include("specials.jl")
include("biomass.jl")
include("finalize.jl")

end
