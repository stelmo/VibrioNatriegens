module Reconstruction 

using DocStringExtensions
using COBREXA

"""
$(TYPEDSIGNATURES)

Entry point function to build the entire curated model.
"""
function build_model()
    model = StandardModel("Vibrio_Natriegens")

    return model
end

end # module