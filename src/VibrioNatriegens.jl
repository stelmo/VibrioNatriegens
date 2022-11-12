"""
VibrioNatriegens

This module exports one function only: `build_model`, which builds the metabolic
model of Vibrio natriegens.
"""
module VibrioNatriegens

using DocStringExtensions

include("ModuleTools.jl")
using .ModuleTools

@inc "Utils"
@inc_dir "metabolism"

using COBREXA.Types

using .Utils
using .Glycolysis_Gluconeogenesis
using .Ala_Asp_Glu

"""
$(TYPEDSIGNATURES)

Entry point function to build the entire curated model.
"""
function build_model()
    model = ObjectModel(id = "Vibrio natriegens")

    # add pathways 
    for rxn in Glycolysis_Gluconeogenesis.rxns
        add_reaction_from_rhea!(model; rxn...)
    end
    
    for rxn in Ala_Asp_Glu.rxns
        add_reaction_from_rhea!(model; rxn...)
    end
    
    return model
end

export build_model

end # module