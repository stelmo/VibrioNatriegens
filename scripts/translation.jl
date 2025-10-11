using DataFrames, CSV, DataFramesMeta



df = DataFrame(CSV.File("data/model/transporters.csv"))
@select!(df, Not(:Substrates))

CSV.write("data/model/transporters.csv", df)
