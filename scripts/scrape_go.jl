using HTTP, JSON
using DataFrames, CSV

function get_descendants(go_id; relationship_types = ["is_a", "part_of", "regulates"])
    relations = join(relationship_types, "%2C")
    url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/$go_id/descendants?relations=$relations" # direct descendents
    resp = HTTP.request("GET", url, ["Accept" => "application/json"])
    resp.status == 200 || error("REST error: $go_id.")
    body = JSON.parse(String(resp.body))
    length(body["results"]) == 1 || error("Multiple results.")
    get(first(body["results"]), "descendants", [])
end

function get_name(go_id)
    url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/$go_id"
    resp = HTTP.request("GET", url, ["Accept" => "application/json"])
    resp.status == 200 || error("REST error: $go_id.")
    body = JSON.parse(String(resp.body))
    length(body["results"]) == 1 || error("Multiple results.")
    get(first(body["results"]), "name", missing)
end

# test
go_id = "GO:0006810" # transport
get_descendants(go_id)
get_name("GO:0015858")

df = DataFrame(GoType = String[], GoSlim = String[], GoID = String[], GoName = String[])
for gotype in ["component", "function", "process"]
    for row in CSV.File(
        joinpath("data", "go_ontology", "goslim_prokaryote_$gotype.txt");
        header = ["GO", "NAME"],
    )
        goslim_id, goslim_name = row.GO, row.NAME
        for desc_id in get_descendants(goslim_id)
            desc_name = get_name(desc_id)
            push!(df, (gotype, goslim_name, desc_id, desc_name))
        end
    end
end

open(joinpath("data", "go_ontology", "goslim_prokaryote_annotations.tab"), "w") do io
    CSV.write(io, df; delim = "\t")
end
