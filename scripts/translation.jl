using DataFrames, CSV, DataFramesMeta

df = DataFrame(CSV.File(joinpath("data", "annotations", "ncbi", "refseq_annotations.tsv")))

trna = @rsubset(df, startswith(:Name, "tRNA"))

ribosomes = @rsubset(df, occursin("ribosomal", :Name))
@subset!(ribosomes, $"Gene Type" .== "protein-coding")

CSV.write("ribosomes.csv", ribosomes)
