using VibrioNatriegens
using AbstractFBCModels
using JSONFBCModels
using SBMLFBCModels

model = VibrioNatriegens.build_model()

VibrioNatriegens.print_reactions(model)
VibrioNatriegens.print_metabolites(model)

m = convert(JSONFBCModels.JSONFBCModel, model)
AbstractFBCModels.save(m, "vnat.json")

m = convert(SBMLFBCModels.SBMLFBCModel, model)
AbstractFBCModels.save(m, "vnat.sbml")

