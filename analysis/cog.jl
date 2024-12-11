using DataFrames, DataFramesMeta, CSV

anno = DataFrame(
    CSV.File(joinpath("data", "annotations", "eggnog", "out.emapper.annotations")),
)
names(anno)
@select!(anno, :query, :eggNOG_OGs)
@rename!(anno, :COG = :eggNOG_OGs)
@rtransform!(anno, :COG = :COG[1:7])
@rsubset!(anno, startswith(:COG, "COG"))

pathways = DataFrame(
    Pathway = String[],
    COG = String[],
    FunctionalCat = String[],
    COG_category = String[],
    Name = String[],
    EC = String[],
)
cog_pathways_base = joinpath("data", "cog", "pathways")
for dir in readdir(joinpath(cog_pathways_base))
    _df = DataFrame(CSV.File(joinpath(cog_pathways_base, dir)))
    @rename!(_df, :FunctionalCat = $"Funtional Cat", :COG_category = $"COG symbol")
    @rtransform!(
        _df,
        :FunctionalCat = :FunctionalCat == false ? "F" : :FunctionalCat,
        :COG_category = :COG_category == false ? "F" : :COG_category
    )
    _df.Pathway = fill(dir[1:end-4], size(_df, 1))
    append!(pathways, _df)
end
@select!(pathways, :Pathway, :COG)

pathways
anno

df = leftjoin(anno, pathways, on = :COG)
dropmissing!(df)

total = @combine(groupby(pathways, :Pathway), :Total = length(unique(:COG)))
leftjoin!(df, total, on = :Pathway)

ddf = @combine(groupby(df, :Pathway), :Completeness = length(unique(:COG)) ./ first(:Total))

@orderby(ddf, :Completeness)
