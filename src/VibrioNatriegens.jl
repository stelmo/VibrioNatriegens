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
using .Citrate_Cycle
using .Pyruvate_Metabolism
using .Pentose_Phosphate_Pathway
using .Pyrimidine_Metabolism
using .Purine_Metabolism
using .Ala_Asp_Glu
using .Cys_Met
using .Val_Leu_Iso_Synthesis
using .Val_Leu_Iso_Degradation
using .Tyr_Metabolism
using .Phe_Metabolism
using .Try_Metabolism
using .Phe_Tyr_Try
using .Gly_Ser_Thr
using .Lysine_Biosynthesis
using .Lysine_Degradation
using .Arigine_Biosynthesis
using .Arg_Pro_Meta
using .Histidine_Meta

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
    
    for rxn in Citrate_Cycle.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Pyruvate_Metabolism.rxns
        add_reaction_from_rhea!(model; rxn...)
    end
    
    for rxn in Pentose_Phosphate_Pathway.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Pyrimidine_Metabolism.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Purine_Metabolism.rxns
        add_reaction_from_rhea!(model; rxn...)
    end
    #flo ^

    for rxn in Ala_Asp_Glu.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Cys_Met.rxns
        add_reaction_from_rhea!(model; rxn...)
    end
    
    for rxn in Val_Leu_Iso_Synthesis.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Val_Leu_Iso_Degradation.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Tyr_Metabolism.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Phe_Metabolism.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Try_Metabolism.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Phe_Tyr_Try.rxns
        add_reaction_from_rhea!(model; rxn...)
    end


    for rxn in Gly_Ser_Thr.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Lysine_Biosynthesis.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Lysine_Degradation.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Arigine_Biosynthesis.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Arg_Pro_Meta.rxns
        add_reaction_from_rhea!(model; rxn...)
    end

    for rxn in Histidine_Meta.rxns
        add_reaction_from_rhea!(model; rxn...)
    end
    return model
end

export build_model

end # module