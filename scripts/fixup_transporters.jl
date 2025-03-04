using CSV, DataFrames, DataFramesMeta

ts = DataFrame(CSV.File(joinpath("data", "annotations", "transaap", "selected_transporters.csv")))

sort(unique(ts.AutoAnnotation))

gdf = groupby(df, [:TransporterFamily, :Substrate])

for gdf in groupby(df, [:TransporterFamily, :Substrate])
    CSV.write(joinpath("data", "annotations", "transaap", "subs", ""))
end

df = DataFrame(Type=String[],Metabolite=String[],CHEBI=String[],Protein=String[],Stoichiometry=String[],Isozyme=String[],)

type_lu = Dict(
    "1.A.11" => "Permease",
    "1.A.8" => "Permease",
    "2.A.1" => "Symport",
    "2.A.13" => "Symport",
    "2.A.14" => "Symport",
    "2.A.15" => "Symport",
    "2.A.20" => "Symport",
    "2.A.21" => "Symport",
    "2.A.23" => "Symport",
    "2.A.25" => "Symport",
    "2.A.26" => "Symport",
    "2.A.27" => "Symport",
    "2.A.28" => "Symport",
    "2.A.3"  => "Symport", # can also be antiport
    "2.A.33" => "Antiport",
    "2.A.34" => "Antiport",
    "2.A.35" => "Antiport",
    "2.A.36" => "Antiport",
    "2.A.37" => "Antiport",
    "2.A.39" => "Symport",
    "2.A.40" => "Symport",
    "2.A.41" => "Symport",
    "2.A.42" => "Symport",
    "2.A.44" => "Symport",
    "2.A.46" => "Symport",
    "2.A.47" => "Symport",
    "2.A.49" => "Antiport",
    "2.A.53" => "Symport",
    "2.A.56" => "Symport",
    "2.A.58" => "Symport",
    "2.A.61" => "Symport", # can also be antiport
    "2.A.62" => "Antiport",
    "2.A.63" => "Antiport",
    "2.A.68" => "Symport",
    "2.A.70" => "Symport",
    "2.A.75" => "Antiport",
    "2.A.76" => "Antiport",
    "2.A.78" => "Antiport",
    "2.A.8" => "Symport",
    "2.A.80" => "Symport",
    "2.A.81" => "Antiport",
    "3.A.1" => "ABC",
    "4.A" => "PTS",
)

@subset(ts, :AutoAnnotation .== "4.A")
