
function add_metabolite_annotations!(model)

    props = [:Chebi, :Kegg, :MolarMass, :Inchi, :Names]
    met_dict = Dict(
        string(row.Chebi) => Dict(
            lowercase(string(p)) =>
                p == :Names ? String.(split(getproperty(row, p), "#")) :
                [string(getproperty(row, p))] for
            p in props if !ismissing(getproperty(row, p))
        ) for row in CSV.File(
            joinpath(
                pkgdir(@__MODULE__),
                "data",
                "annotations",
                "metabolite_annotations.csv",
            ),
        )
    )

    for _mid in A.metabolites(model)

        mid = join(first(split(_mid, "_"))) # this is the core mid
        annos = get!(met_dict, mid, Dict{String,Vector{String}}())

        if RheaReactions.is_cached("metabolites", mid) # glycogen causes this :(
            m = get_metabolite(mid) # overwrite annotations with data from Rhea
            isnothing(m.inchikey) || (annos["inchikey"] = [m.inchikey])
            isnothing(m.smiles) || (annos["smiles"] = [m.smiles])
            isnothing(m.inchi) || (annos["inchi"] = [m.inchi])
        end

        # chebi is the default name space, ensure it is added
        if haskey(annos, "chebi")
            if mid ∉ annos["chebi"]
                push!(annos["chebi"], mid)
            end
        else
            annos["chebi"] = [mid]
        end
        annos["SBO"] = ["SBO_0000299"]

        model.metabolites[_mid].annotations = annos
    end

    model.metabolites["glycogen"].annotations = Dict(
        "names" => ["glycogen"],
        "kegg" => ["C00182"],
        "molarmass" => ["162.1406"],
        "SBO" => ["SBO_0000299"],
        "chebi" => ["28087"],
    )
end
