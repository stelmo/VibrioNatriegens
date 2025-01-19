using VibrioNatriegens
using AbstractFBCModels

model = VibrioNatriegens.build_model()

VibrioNatriegens.print_reactions(model)
VibrioNatriegens.print_metabolites(model)
