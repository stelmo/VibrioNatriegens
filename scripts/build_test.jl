using VibrioNatriegens
using AbstractFBCModels
using JSONFBCModels
using SBMLFBCModels
import AbstractFBCModels as A
using RheaReactions, DataFrames, DataFramesMeta, CSV

model = VibrioNatriegens.build_model()

VibrioNatriegens.print_reactions(model)
VibrioNatriegens.print_metabolites(model)

m = convert(JSONFBCModels.JSONFBCModel, model)
AbstractFBCModels.save(m,"vnat.json")

m = convert(SBMLFBCModels.SBMLFBCModel, model)
AbstractFBCModels.save(m,"vnat.sbml")

