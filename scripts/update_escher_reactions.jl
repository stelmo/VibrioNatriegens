using DataFrames, DataFramesMeta, CSV, JSON, VibrioNatriegens

escher_maps = [
    "amino_acids.json"
    "carbohydrate_metabolism.json"
    "cofactors_etc.json"
    "energy_metabolism.json"
    "lipids.json"
    "nucleotides.json"
]

m = JSON.parsefile(joinpath("maps", escher_maps[1]))
m[2]["reactions"]["1"]


model = VibrioNatriegens.build_model()
