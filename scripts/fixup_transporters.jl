using CSV, DataFrames, DataFramesMeta

ts = DataFrame(
    CSV.File(joinpath("data", "annotations", "transaap", "selected_transporters.csv")),
)

sort(unique(ts.AutoAnnotation))

gdf = groupby(df, [:TransporterFamily, :Substrate])

for gdf in groupby(df, [:TransporterFamily, :Substrate])
    CSV.write(joinpath("data", "annotations", "transaap", "subs", ""))
end

df = DataFrame(
    Type = String[],
    Metabolite = String[],
    CHEBI = String[],
    Protein = String[],
    Stoichiometry = String[],
    Isozyme = String[],
)



gs = [
    "WP_020333513.1"
    "WP_014233759.1"
    "WP_020336058.1"
    "WP_020336067.1"
]
f = open("temp.txt", "w")
for (k, kk) in [
    ("leucine", 57427),
    ("isoleucine", 58045),
    ("valine", 57762),
    ("cysteine", 35235),
    ("alanine", 57972),
    ("serine", 33384),
    ("phenylalanine", 58095),
    ("tyrosine", 58315),
    ("glutamate", 29985),
    ("glycine", 57305),
    ("asparate", 29991),
    ("lysine", 32551),
    ("methionine", 57844),
    ("tryptophan", 57912),
    ("histidine", 57595),
    ("asparagine", 58048),
    ("arginine", 32682),
    ("threonine", 57926),
    ("glutamine", 58359),
    ("proline", 60039),
]
    for (i, g) in enumerate(gs)
        c = cs[i]
        write(f, "Symport,$k/H,CHEBI:$kk/CHEBI:15378,$g,1,$c\n")
    end
end
close(f)

