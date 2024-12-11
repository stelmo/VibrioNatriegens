using CSV, DataFrames, DataFramesMeta, ReadableRegex


seed = DataFrame(CSV.File(joinpath("data", "seed", "reactions.tsv")))
@select!(seed, :id, :aliases)
reg = @compile exactly(1, "KEGG: R") * exactly(5, DIGIT)
findkegg(entry) = begin
    x = [string(last(split(m.match, " "))) for m in eachmatch(reg, entry)]
    isempty(x) ? missing : x
end
@rtransform!(seed, :KEGG = findkegg.(:aliases))
dropmissing!(seed)
seed = flatten(seed, :KEGG)
@rename!(seed, :Reaction = :id)
@select!(seed, :Reaction, :KEGG)

df = DataFrame(CSV.File(joinpath("data", "annotations", "grr", "updated_grrs.txt")))
@rtransform!(df, :Reaction = :Reaction[3:end-2])

mdf = leftjoin(df, seed, on = :Reaction)
dropmissing!(mdf)
mdf = @select mdf Not(:Reaction)
@rename!(mdf, :Reaction = :KEGG)

# rename proteins to new

old = DataFrame(CSV.File(joinpath("data", "annotations", "ncbi", "ncbi_annotations.tsv")))
@rename!(old, :LocusTag = $"Locus tag")
@select!(old, :Accession, :Begin, :End, :Chromosome, :Orientation, :Name, :LocusTag)

new = DataFrame(CSV.File(joinpath("data", "annotations", "ncbi", "refseq_annotations.tsv")))
@rename!(new, :ProteinAccession = $"Protein accession")
@select!(new, :Accession, :Begin, :End, :Chromosome, :Orientation, :Name, :ProteinAccession)
@rtransform!(new, :Accession = :Accession[4:end])

df = innerjoin(
    new,
    old,
    on = [:Accession, :Begin, :Chromosome, :Orientation],
    makeunique = true,
)
@select!(df, :ProteinAccession, :LocusTag)
dropmissing!(df)
@rtransform!(df, :LocusTag = replace(:LocusTag, "_" => ""))
@rename!(df, :Gene = :LocusTag)

leftjoin!(mdf, df, on = :Gene)
dropmissing!(mdf)

CSV.write(joinpath("data", "annotations", "grr", "mapped_grrs.csv"), mdf)

