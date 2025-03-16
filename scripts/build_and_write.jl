using VibrioNatriegens
using AbstractFBCModels
using JSONFBCModels
using SBMLFBCModels

model = VibrioNatriegens.build_model()

VibrioNatriegens.print_reactions(model)
VibrioNatriegens.print_metabolites(model)

m = convert(JSONFBCModels.JSONFBCModel, model)
AbstractFBCModels.save(m, "vibrio_natriegens.json")

m = convert(SBMLFBCModels.SBMLFBCModel, model)
AbstractFBCModels.save(m, "vibrio_natriegens.sbml")

