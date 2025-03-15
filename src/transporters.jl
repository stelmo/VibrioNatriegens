
"""
$(TYPEDSIGNATURES)

Add exchange reactions for basic metabolites
"""
function add_exchanges!(model)

    df = DataFrame(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "model", "exchange_metabolites.csv"),
            types = [String, String],
        ),
    )

    # by default, exchanges can only export metabolites
    for mid in df.CHEBI

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

    # all extracellular move with diffusion into periplasm
    df = DataFrame(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "model", "exchange_metabolites.csv"),
            types = [String, String],
        ),
    )

    # bidirectional by default
    for mid in df.CHEBI

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

    df = DataFrame(
        CSV.File(joinpath(pkgdir(@__MODULE__), "data", "model", "transporters.csv")),
    )

    gs = String[]
    ms = String[]

    # abc transporters
    abcs = @subset(df, :Type .== "ABC")
    for g in groupby(abcs, [:CHEBI, :Isozyme])
        mid = first(g.CHEBI)
        if mid in A.metabolites(model)
            push!(ms, mid)
            iso = string.(g.Protein)
            ss = parse.(Float64, string.(g.Stoichiometry))
            append!(gs, iso)
            add_abc!(model, mid, iso, ss)
        else
            @warn "$mid not in model (ABC)"
        end
    end

    # PTS transporters
    pts = @subset(df, :Type .== "PTS")
    for g in groupby(pts, [:CHEBI, :Isozyme])
        mid = first(g.CHEBI)
        if mid in A.metabolites(model)
            push!(ms, mid)
            iso = string.(g.Protein)
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_pts!(model, mid, iso, ss)
        else
            @warn "$mid not in model (PTS)"
        end
    end

    # symport
    symport = @subset(df, :Type .== "Symport")
    for g in groupby(symport, [:CHEBI, :Isozyme])
        # println(g)
        mid1, mid2 = sort(split(first(g.CHEBI), "/")) # to make rid unique
        if mid1 in A.metabolites(model) && mid2 in A.metabolites(model)
            push!(ms, mid1)
            push!(ms, mid2)
            iso = string.(g.Protein)
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_symport!(model, mid1, mid2, iso, ss)
        else
            @warn "$mid1 or $mid2 not in model (symport)"
        end
    end

    # antiport
    antiport = @subset(df, :Type .== "Antiport")
    for g in groupby(antiport, [:CHEBI, :Isozyme])
        mid1, mid2 = sort(split(first(g.CHEBI), "/")) # to make rid unique
        if mid1 in A.metabolites(model) && mid2 in A.metabolites(model)
            push!(ms, mid1)
            push!(ms, mid2)
            iso = string.(g.Protein)
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_antiport!(model, mid1, mid2, iso, ss)
        else
            @warn "$mid1 or $mid2 not in model (antiport)"
        end
    end

    # permease (the default as well)
    permease = @subset(df, :Type .== "Permease")
    for g in groupby(permease, [:CHEBI, :Isozyme])
        mid = first(g.CHEBI)
        if mid in A.metabolites(model)
            push!(ms, mid)
            iso = string.(g.Protein)
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_permease!(model, mid, iso, ss)
        else
            @warn "$mid not in model (permease)"
        end
    end

    # add default permeases - only reactions that did not get another transporter
    all_exchange_metabolites =
        DataFrame(
            CSV.File(
                joinpath(pkgdir(@__MODULE__), "data", "model", "exchange_metabolites.csv"),
            ),
        ).CHEBI
    # note: Pi, Na, and H will not get a permease here, due to them being involved in the other porters
    missing_transporters = setdiff(string.(all_exchange_metabolites), unique(ms))
    for mid in missing_transporters
        if mid in A.metabolites(model)
            add_permease!(model, mid, nothing, [1.0])
        end
    end

    add_genes!(model, gs)
    # no need to add metabolites, because they should all already be in the model
    @assert all(in.(ms, Ref(A.metabolites(model))))
end

function add_abc!(model, mid, iso, ss)
    rid = "ABC_$mid"
    isoz =
        isnothing(iso) ? nothing :
        X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        isnothing(isoz) || push!(model.reactions[rid].gene_association, isoz)
    else
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
            gene_association = isnothing(isoz) ? nothing : [isoz],
            annotations = Dict("SBO" => ["SBO_0000284"]),
        )
    end
end

function add_pts!(model, mid, iso, ss)
    lu_phospho = Dict(
        "506227" => "57513", # n-acetyl-glucosamine -> N-acetyl-D-glucosamine 6-phosphate
        "15903" => "58247", # glucose -> glucose 6 phosphate
        "17992" => "57723", # sucrose -> sucrose 6 phosphate
        "16899" => "61381", # mannitol -> D-mannitol 1-phosphate
        "28645" => "57634", # β-D-fructose -> beta-D-fructose 6-phosphate
    )

    rid = "PTS_$mid"
    isoz =
        isnothing(iso) ? nothing :
        X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        isnothing(isoz) || push!(model.reactions[rid].gene_association, isoz)
    else
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
            gene_association = isnothing(isoz) ? nothing : [isoz],
            annotations = Dict("SBO" => ["SBO_0000284"]),
        )
    end
end

function add_symport!(model, mid1, mid2, iso, ss, dir=:bidir)
    if dir == :forward
        lb = 0.0
        ub = 1000.0    
    elseif dir == :reverse
        lb = -1000.0
        ub = 0.0    
    else
        lb = -1000.0
        ub = 1000.0    
    end

    rid = "SYM_$(mid1)_$mid2"
    isoz =
        isnothing(iso) ? nothing :
        X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        isnothing(isoz) || push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Symport $(A.metabolite_name(model, String(mid1)))::$(A.metabolite_name(model, String(mid2)))",
            stoichiometry = Dict(
                mid1 * "_p" => -1.0,
                mid1 => 1.0,
                mid2 * "_p" => -1.0,
                mid2 => 1.0,
            ),
            objective_coefficient = 0.0,
            lower_bound = lb,
            upper_bound = ub,
            gene_association = isnothing(isoz) ? nothing : [isoz],
            annotations = Dict("SBO" => ["SBO_0000284"]),
        )
    end
end

function add_antiport!(model, mid1, mid2, iso, ss, dir=:bidir)

    if dir == :forward
        lb = 0.0
        ub = 1000.0    
    elseif dir == :reverse
        lb = -1000.0
        ub = 0.0    
    else
        lb = -1000.0
        ub = 1000.0    
    end

    rid = "ANTI_$(mid1)_$mid2"
    isoz =
        isnothing(iso) ? nothing :
        X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        isnothing(isoz) || push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Antiport $(A.metabolite_name(model, String(mid1)))::$(A.metabolite_name(model, String(mid2)))",
            stoichiometry = Dict(
                mid1 * "_p" => 1.0,
                mid1 => -1.0,
                mid2 * "_p" => -1.0,
                mid2 => 1.0,
            ),
            objective_coefficient = 0.0,
            lower_bound = lb,
            upper_bound = ub,
            gene_association = isnothing(isoz) ? nothing : [isoz],
            annotations = Dict("SBO" => ["SBO_0000284"]),
        )
    end
end

function add_permease!(model, mid, iso, ss)
    rid = "PERM_$mid"
    isoz =
        isnothing(iso) ? nothing :
        X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        isnothing(isoz) || push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Permease $(A.metabolite_name(model, String(mid)))",
            stoichiometry = Dict(mid * "_p" => -1.0, mid => 1.0),
            objective_coefficient = 0.0,
            lower_bound = -1000,
            upper_bound = 1000,
            gene_association = isnothing(isoz) ? nothing : [isoz],
            annotations = Dict("SBO" => ["SBO_0000284"]),
        )
    end
end
