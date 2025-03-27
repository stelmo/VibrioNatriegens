using CSV, DataFrames, DataFramesMeta

df = DataFrame(CSV.File(joinpath("data", "annotations", "kegg", "ko.txt"), header=["Gene","KO","Desc"]))
