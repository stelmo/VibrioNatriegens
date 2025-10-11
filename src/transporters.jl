
"""
$(TYPEDSIGNATURES)

Add exchange reactions for basic metabolites
"""
function add_exchanges!(model)
    # pkgdir(@__MODULE__)

    # by default, exchanges can only export metabolites
    for row in CSV.File(
        joinpath(pkgdir(@__MODULE__), "data", "model", "exchanges.csv"),
        types = [String, String],
    )
        mid = row.Chebi
        nm = row.Name

        nm = A.metabolite_name(model, mid)

        model.metabolites[mid*"_e"] = deepcopy(model.metabolites[mid])

        model.metabolites[mid*"_e"].compartment = "Extracellular"

        model.reactions["EX_"*mid] = Reaction(;
            name = "Exchange $nm",
            stoichiometry = Dict(mid * "_e" => -1),
            lower_bound = 0,
            upper_bound = 1000,
            annotations = Dict("SBO" => ["SBO_0000284", "SBO_0000627"]),
        )
    end
end

function add_periplasm_transporters!(model)

    # all extracellular metabolites move with diffusion into periplasm
    # pkgdir(@__MODULE__), 

    # bidirectional by default
    for row in CSV.File(
        joinpath(pkgdir(@__MODULE__), "data", "model", "exchanges.csv"),
        types = [String, String],
    )
        mid = row.Chebi

        nm = A.metabolite_name(model, mid)

        model.metabolites[mid*"_p"] = deepcopy(model.metabolites[mid])

        model.metabolites[mid*"_p"].compartment = "Periplasm"

        model.reactions["DF_"*mid] = Reaction(;
            name = "Diffusion $nm",
            stoichiometry = Dict(mid * "_e" => -1, mid * "_p" => 1),
            lower_bound = -1000.0,
            upper_bound = 1000.0,
            annotations = Dict("SBO" => ["SBO_0000284"]),
        )
    end

end

function add_membrane_transporters!(model)

    # pkgdir(@__MODULE__), 
    acronym = Dict(
        "Permease" => "PERM",
        "ABC" => "ABC",
        "Symport" => "SYM",
        "Antiport" => "ANTI",
        "PTS" => "PTS",
    )
    gs = String[]
    ms = String[]

    for row in CSV.File(joinpath(pkgdir(@__MODULE__), "data", "model", "transporters.csv"))
        type = row.Type
        protein = String(row.Protein)
        push!(gs, protein)
        stoich = row.Stoichiometry
        iso = String(row.Isozyme)

        mids = String.(sort(split(row.Chebi, "/"))) # order metabolites
        @assert all(haskey(model.metabolites, mid) for mid in mids) "Metabolite(s) $mids not in model"
        urid = acronym[type] * "_" * join(mids, "_") # unique id for reaction given type and metabolites transported

        if haskey(model.reactions, urid) # extend the gene reaction rule

            isos = model.reactions[urid].gene_association
            if haskey(isos, iso) # add to complex
                isos[iso].gene_product_stoichiometry[protein] = stoich
            else # new isozyme
                isos[iso] = Isozyme(Dict(protein => stoich), nothing, nothing)
            end

        else # add a new transporter reaction

            if type == "Permease"
                add_permease!(model, urid, first(mids))
            elseif type == "ABC"
                add_abc!(model, urid, first(mids))
            elseif type == "Symport"
                add_symport!(model, urid, mids...)
            elseif type == "Antiport"
                add_antiport!(model, urid, mids...)
            elseif type == "PTS"
                add_pts!(model, urid, first(mids))
            else
                @warn "Unknown transporter type encountered"
            end

            # add a new gene reaction rule to transport reaction
            model.reactions[urid].gene_association =
                Dict(iso => Isozyme(Dict(protein => stoich), nothing, nothing))
        end

        # push to list of metabolites that have a transporter added
        append!(ms, mids)
    end

    # add default permeases - only reactions that did not get another transporter
    all_exchange_metabolites = CSV.File(
        joinpath(pkgdir(@__MODULE__), "data", "model", "exchanges.csv"),
        types = [String, String],
    ).Chebi
    # note: Pi, Na, and H will not get a permease here, due to them being involved in the other porters

    missing_transporters = setdiff(string.(all_exchange_metabolites), unique(ms))

    for mid in missing_transporters
        if mid in A.metabolites(model)
            add_permease!(model, "PERM_$mid", mid) # these reactions won't have a gene reaction rule, since they are fallbacks
        end
    end

    add_genes!(model, gs)
end

function add_abc!(model, rid, mid)
    st = Dict(
        "30616" => -1, # atp
        "15377" => -1, # water
        "43474" => 1, # pi
        "456216" => 1, # adp
        "15378" => 1,  # h+ 
        mid * "_p" => -1.0,
    )
    st[mid] = get(st, mid, 0) + 1.0 # handle the case when phosphate is transported
    model.reactions[rid] = Reaction(
        name = "Transport $(A.metabolite_name(model, String(mid))) ABC",
        stoichiometry = st,
        objective_coefficient = 0.0,
        lower_bound = 0,
        upper_bound = 1000,
        gene_association = nothing,
        annotations = Dict("SBO" => ["SBO_0000284"]),
    )
end

function add_pts!(model, rid, mid)
    lu_phospho = Dict(
        "506227" => "57513", # n-acetyl-glucosamine -> N-acetyl-D-glucosamine 6-phosphate
        "15903" => "58247", # glucose -> glucose 6 phosphate
        "17992" => "57723", # sucrose -> sucrose 6 phosphate
        "16899" => "61381", # mannitol -> D-mannitol 1-phosphate
        "28645" => "57634", # β-D-fructose -> beta-D-fructose 6-phosphate
    )
    model.reactions[rid] = Reaction(
        name = "Transport $(A.metabolite_name(model, String(mid))) PTS",
        stoichiometry = Dict(
            "58702" => -1.0, # pep
            "15361" => 1.0, # pyr
            mid * "_p" => -1.0,
            lu_phospho[mid] => 1.0, # cytosol phospho metabolite
        ),
        objective_coefficient = 0.0,
        lower_bound = 0,
        upper_bound = 1000,
        gene_association = nothing,
        annotations = Dict("SBO" => ["SBO_0000284"]),
    )
end

function add_symport!(model, rid, mid1, mid2)
    model.reactions[rid] = Reaction(
        name = "Symport $(A.metabolite_name(model, String(mid1)))::$(A.metabolite_name(model, String(mid2)))",
        stoichiometry = Dict(
            mid1 * "_p" => -1.0,
            mid1 => 1.0,
            mid2 * "_p" => -1.0,
            mid2 => 1.0,
        ),
        objective_coefficient = 0.0,
        lower_bound = -1000,
        upper_bound = 1000,
        gene_association = nothing,
        annotations = Dict("SBO" => ["SBO_0000284"]),
    )
end

function add_antiport!(model, rid, mid1, mid2)
    model.reactions[rid] = Reaction(
        name = "Antiport $(A.metabolite_name(model, String(mid1)))::$(A.metabolite_name(model, String(mid2)))",
        stoichiometry = Dict(
            mid1 * "_p" => 1.0,
            mid1 => -1.0,
            mid2 * "_p" => -1.0,
            mid2 => 1.0,
        ),
        objective_coefficient = 0.0,
        lower_bound = -1000,
        upper_bound = 1000,
        gene_association = nothing,
        annotations = Dict("SBO" => ["SBO_0000284"]),
    )
end

function add_permease!(model, rid, mid)
    model.reactions[rid] = Reaction(
        name = "Permease $(A.metabolite_name(model, String(mid)))",
        stoichiometry = Dict(mid * "_p" => -1.0, mid => 1.0),
        objective_coefficient = 0.0,
        lower_bound = -1000,
        upper_bound = 1000,
        gene_association = nothing,
        annotations = Dict("SBO" => ["SBO_0000284"]),
    )
end
