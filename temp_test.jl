using VibrioNatriegens, COBREXA, Gurobi


model = VibrioNatriegens.build_model(;pathway_maps=["map00010"]
)
