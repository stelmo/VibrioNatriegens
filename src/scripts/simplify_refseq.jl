using CSV, DataFrames, DataFramesMeta, ReadableRegex

new = DataFrame(CSV.File(joinpath("data", "annotations", "ncbi", "refseq_annotations.tsv")))
@select!(new, $"Protein accession", :Name, $"Locus tag")


CSV.write(joinpath("refseq_annotations_simplified.csv"), new)

