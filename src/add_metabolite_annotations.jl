function add_metabolite_annotations!(model)

    names_df = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "chebi",
                "names.tsv",
            ),
        ),
    )
    @subset!(names_df, :LANGUAGE .== "en")
    names_dict = Dict(
        "CHEBI:"*string(first(gdf.COMPOUND_ID)) =>
            unique(gdf.NAME)
         for gdf in groupby(names_df, :COMPOUND_ID)
    )

    inchi_df = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "chebi",
                "chebiId_inchi.tsv",
            ),
        ),
    )
    inchi_dict = Dict(zip("CHEBI:".*string.(inchi_df.CHEBI_ID), inchi_df.InChI))

    chemical_df = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "chebi",
                "chemical_data.tsv",
            ),
        ),
    )
    @subset!(chemical_df, :SOURCE .== "ChEBI", :TYPE .== "MASS")
    chemical_dict = Dict(
        "CHEBI:"*string(first(gdf.COMPOUND_ID)) => parse(Float64, first(gdf.CHEMICAL_DATA))
         for gdf in groupby(chemical_df, :COMPOUND_ID)
    )

    accession_df = DataFrame(
        CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "chebi",
                "database_accession.tsv",
            ),
        ),
    )
    kegg_compound = Dict(
        "CHEBI:"*string(first(gdf.COMPOUND_ID)) => gdf.ACCESSION_NUMBER
        for gdf in groupby(@subset(accession_df, :TYPE .== "KEGG COMPOUND accession"), :COMPOUND_ID)
    )
 
    for mid in A.metabolites(model)
        m = model.metabolites[mid] # name, charge, formula, compartment already added
        m.molarmass = get(chemical_dict, mid, nothing)
        m.annotations = Dict(
            "Names" => get(names_dict, mid, ["",]),
            "InChI" => [get(inchi_dict, mid, ""),],
            "KEGG" => get(kegg_compound, mid, ["",]),
        )
    end


    model.metabolites["glycogen"].molarmass = 162.1406
    model.metabolites["glycogen"].annotations = Dict(
        "Names" => ["glycogen",],
        "InChI" => ["",],
        "KEGG" => ["C00182",],
    )

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
