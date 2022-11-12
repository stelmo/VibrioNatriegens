module Utils

using ..ModuleTools

using DocStringExtensions
using HTTP
using COBREXA
using COBREXA.Types
using COBREXA.Reconstruction
import RheaReactions as rr

"""
$(TYPEDSIGNATURES)

Add a reaction with Rhea ID `rhea_id` to `model`. Throw an error if already
present. Here `grr` refers to the uniprot IDs. This is corrected to the gene
names in a later processing step.
"""
function add_reaction_from_rhea!(
    model;
    rhea_id::Integer,
    name = "",
    isozymes=nothing,
    subsystem=nothing,
)
    rxn = rr.get_reaction(rhea_id)
    rxn.accession in reactions(model) && throw(error("Reaction with ID $rhea_id already in model."))
    !rxn.isbalanced && @warn("Reaction not balanced.")
    coeff_mets = rr.get_reaction_metabolites(rhea_id)
    metabolite_dict = Dict(m.accession => s for (s, m) in coeff_mets)

    # add gene to model
    if !isnothing(isozymes)
        for isozyme in isozymes 
            for (gid, _) in isozyme
                gid ∉ genes(model) && add_gene!(model, Gene(gid))
            end
        end    
    end
    
    # add reaction to model 
    isos = isnothing(isozymes) ?  nothing : [Isozyme(stoichiometry = Dict(k => v for (k, v) in iso)) for iso in isozymes]
    add_reaction!(
        model, 
        Reaction(
            rxn.accession;
            name,
            metabolites=metabolite_dict,
            gene_associations=isos,
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

@export_locals()

end # module