module Reconstruction 

using DocStringExtensions
using HTTP
using COBREXA
import KEGGModuleParser as kp
import BiGGReactions as br
import RheaReactions as rr
import MetaNetXReactions as mnr

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

Add a reaction with Rhea ID `rhea_id` to `model`. Throw an error if already
present. Here `grr` refers to the uniprot IDs. This is corrected to the gene
names in a later processing step.
"""
function add_reaction_from_rhea!(
    model, 
    rhea_id::Integer; 
    name=nothing, 
    grrs=nothing,
    subsystem=nothing,
)
    rhea_id in reactions(model) && throw(error("Reaction with ID $rhea_id already in model."))

    rxn = rr.get_reaction(rhea_id)
    !rxn.isbalanced && @warn("Reaction not balanced.")
    coeff_mets = rr.get_reaction_metabolites(rhea_id)
    metabolite_dict = Dict(m.accession => s for (s, m) in coeff_mets)

    # add gene to model
    for grr in grrs 
        for g in grr 
            g ∉ genes(model) && add_gene!(model, Gene(g))
        end
    end
    
    # add reaction to model 
    add_reaction!(
        model, 
        Reaction(
            rxn.accession;
            name,
            metabolites=metabolite_dict,
            grr=grrs,
            subsystem,
        ),
    )

    # add metabolites to model
    for (_, met) in coeff_mets 
        if met.accession ∉ metabolites(model)
            add_metabolite!(
                model, 
                Metabolite(
                    met.accession;
                    name = met.name,
                    formula = met.formula,
                    charge = met.charge,
                ),
            )
        end
    end

    nothing
end


"""
$(TYPEDSIGNATURES)

Get all the EC numbers associated with KEGG module `mnum`.
"""
get_ecs(mnum) = last.(split.(kp.get_module_ECs(mnum),":"))

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