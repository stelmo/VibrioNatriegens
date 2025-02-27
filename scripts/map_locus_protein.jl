using FASTX, CSV, DataFrames

df = DataFrame(LocusTag=String[], ProteinID=String[])

FASTAReader(open(joinpath("data", "genome", "proteome.fasta"))) do reader
    for record in reader
        lt = string(first(split(string(last(split(description(record), "locus_tag="))), "]")))
        pid = string(first(split(string(last(split(description(record), "protein_id="))), "]")))
        length(pid) <= 20 && push!(df, (lt, pid))
    end
end

CSV.write(joinpath("data", "genome", "locustag_proteinid.csv"), df)
