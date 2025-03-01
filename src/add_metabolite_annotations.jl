function add_metabolite_annotations!(model)

    met_df = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "chebi",
                "metabolite_annotations.csv",
            ),
        ),
    )

    met_dict = Dict(
        string(first(gdf.MetaboliteID)) => Dict(
            k => String(string(v)) for
            (k, v) in zip(labels(gdf), gdf[1, :]) if !ismissing(v)
        ) for gdf in groupby(met_df, :MetaboliteID)
    )

    for mid in A.metabolites(model)
        _d = get(met_dict, mid, Dict())
        m = model.metabolites[mid] # name, charge, formula, compartment already added
        m.molarmass = parse(Float64, get(_d, "MolarMass", "0"))
        m.annotations = Dict(
            "Names" => string.(split(get(_d, "Names", ""), "#")),
            "InChI" => [get(_d, "InChI", ""),],
            "KEGG" => string.(split(get(_d, "KeGG", ""),"#")),
        )
    end

    model.metabolites["glycogen"].molarmass = 162.1406
    model.metabolites["glycogen"].annotations =
        Dict("Names" => ["glycogen"], "InChI" => [""], "KEGG" => ["C00182"])

end

function fix_noncytosolic_metabolite_annotations!(model)

    for mid in A.metabolites(model)
        if occursin("_p", mid) || occursin("_e", mid)
            m = model.metabolites[mid]
            fully_annotated = model.metabolites[first(split(mid, "_"))] # name, charge, formula, compartment already added
            m.molarmass = fully_annotated.molarmass
            m.annotations = fully_annotated.annotations
            m.name = fully_annotated.name
            m.charge = fully_annotated.charge
            m.formula = fully_annotated.formula
        end
    end
end
