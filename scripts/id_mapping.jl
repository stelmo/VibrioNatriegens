using CSV, DataFrames, DataFramesMeta
using VibrioNatriegens, JSON

model = VibrioNatriegens.build_model()
mids = collect(keys(model.metabolites))

df = DataFrame(CSV.File(joinpath("data", "experiments", "coppens_2023_biolog.csv")))

j1 = JSON.parsefile("cpds_1.json")
j2 = JSON.parsefile("cpds_2.json")

lus = []
for x in df.ID
    xmids = haskey(j1, x) ? j1[x]["xrefs"] : j2[x]["xrefs"]
    idxs = findall(in.(xmids, Ref(mids)))
    push!(lus, isempty(idxs) ? missing : first(xmids[idxs]))
end
lus

df.Chebi = lus

CSV.write(joinpath("data", "experiments", "coppens_2023_biolog.csv"), df)



