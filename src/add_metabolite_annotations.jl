function add_metabolite_annotations!(model)

    met_df = DataFrame(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "chebi", "metabolite_annotations.csv"),
        ),
    )
    renamer(k) = begin
        if k == "kegg"
            return "kegg.compound"
        else
            return k
        end
    end
    met_dict = Dict(
        string(first(gdf.chebi)) => Dict(
            renamer(k) => [String(string(v))] for
            (k, v) in zip(labels(gdf), gdf[1, :]) if !ismissing(v)
        ) for gdf in groupby(met_df, :chebi)
    )

    for (k, v) in met_dict
        if haskey(v, "names")
            met_dict[k]["names"] = string.(split(first(v["names"]), "#"))
        end
    end

    for _mid in A.metabolites(model)

        mid = join(first(split(_mid, "_"))) # this is the basic mid
        mm = get!(met_dict, mid, Dict{String,String}())

        if occursin("_p", _mid) || occursin("_e", _mid)
            mm = get!(met_dict, _mid, mm)
        end

        if RheaReactions.is_cached("metabolites", mid) # glycogen causes this :(
            m = get_metabolite(mid)
            isnothing(m.inchikey) || (mm["inchikey"] = [m.inchikey])
            isnothing(m.smiles) || (mm["smiles"] = [m.smiles])
            isnothing(m.inchi) || (mm["inchi"] = [m.inchi])
        end

        # chebi is the default name space, ensure it is added
        if haskey(mm, "chebi")
            if mid ∉ mm["chebi"]
                push!(mm["chebi"], mid)
            end
        else
            mm["chebi"] = [mid]
        end
        mm["SBO"] = ["SBO_0000299"]
    end

    for mid in A.metabolites(model)
        m = model.metabolites[mid] # name, charge, formula, compartment already added
        m.annotations = get(met_dict, mid, Dict{String,Vector{String}}())
    end

    model.metabolites["glycogen"].annotations = Dict(
        "names" => ["glycogen"],
        "kegg" => ["C00182"],
        "molarmass" => ["162.1406"],
        "SBO" => ["SBO_0000299"],
        "chebi" => ["28087"],
    )
end
