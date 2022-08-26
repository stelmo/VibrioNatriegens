module Reconstruction 

using DocStringExtensions
using HTTP
using COBREXA

"""
$(TYPEDSIGNATURES)

Entry point function to build the entire curated model.
"""
function build_model()
    model = StandardModel("Vibrio_Natriegens")

    return model
end

"""
$(TYPEDSIGNATURES)

Downloads reference models and stores them in `savedir`.
"""
function download_models(;savedir=joinpath("data", "models"))
    iml1515_path = joinpath("data", "models", "iml1515.json") 
    isfile(iml1515_path) || download("http://bigg.ucsd.edu/static/models/iML1515.json", iml1515_path)

    return nothing
end

end # module