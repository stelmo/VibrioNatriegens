
"""
$(TYPEDEF)

Reaction field adapted to enzyme constrained models. 
`gene_association` contains `COBREXA.Isozyme`s.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Reaction
    name::A.Maybe{String} = nothing
    lower_bound::Float64 = -Inf
    upper_bound::Float64 = Inf
    stoichiometry::Dict{String,Float64} = Dict()
    objective_coefficient::Float64 = 0.0
    gene_association::A.Maybe{Vector{X.Isozyme}} = nothing
    dg::A.Maybe{Float64} = nothing
    annotations::A.Annotations = A.Annotations()
    notes::A.Notes = A.Notes()
    transporter::Bool = false
end

Base.show(io::Base.IO, ::MIME"text/plain", x::Reaction) = A.pretty_print_kwdef(io, x)

"""
$(TYPEDEF)

Normal metabolite.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Metabolite
    name::A.Maybe{String} = nothing
    compartment::A.Maybe{String} = nothing
    formula::A.Maybe{A.MetaboliteFormula} = nothing
    charge::A.Maybe{Int} = nothing
    molarmass::A.Maybe{Float64} = nothing
    balance::Float64 = 0.0
    annotations::A.Annotations = A.Annotations()
    notes::A.Notes = A.Notes()
end

Base.show(io::Base.IO, ::MIME"text/plain", x::Metabolite) = A.pretty_print_kwdef(io, x)

"""
$(TYPEDEF)

Gene with an additional molar mass field.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Gene
    name::A.Maybe{String} = nothing
    symbol::A.Maybe{String} = nothing
    molarmass::A.Maybe{Float64} = nothing
    annotations::A.Annotations = A.Annotations()
    notes::A.Notes = A.Notes()
end

Base.show(io::Base.IO, ::MIME"text/plain", x::Gene) = A.pretty_print_kwdef(io, x)

"""
$(TYPEDEF)

Custom model type to represent _Vibrio natriegens_.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Model <: A.AbstractFBCModel
    reactions::Dict{String,Reaction} = Dict()
    metabolites::Dict{String,Metabolite} = Dict()
    genes::Dict{String,Gene} = Dict()
    couplings::Dict{String,A.CanonicalModel.Coupling} = Dict()
end

Base.show(io::Base.IO, ::MIME"text/plain", model::Model) = println("""
    Vibrio natriegens
    $(A.n_reactions(model)) reactions
    $(A.n_metabolites(model)) metabolites
    $(A.n_genes(model)) genes
""")

A.reactions(m::Model) = sort(collect(keys(m.reactions)))
A.metabolites(m::Model) = sort(collect(keys(m.metabolites)))
A.genes(m::Model) = sort(collect(keys(m.genes)))
A.couplings(m::Model) = sort(collect(keys(m.couplings)))
A.n_reactions(m::Model) = length(m.reactions)
A.n_metabolites(m::Model) = length(m.metabolites)
A.n_genes(m::Model) = length(m.genes)
A.n_couplings(m::Model) = length(m.couplings)
A.reaction_name(m::Model, id::String) = m.reactions[id].name
A.metabolite_name(m::Model, id::String) = m.metabolites[id].name
A.gene_name(m::Model, id::String) = m.genes[id].name
A.coupling_name(m::Model, id::String) = m.couplings[id].name
A.reaction_annotations(m::Model, id::String) = m.reactions[id].annotations
A.metabolite_annotations(m::Model, id::String) = m.metabolites[id].annotations
A.gene_annotations(m::Model, id::String) = m.genes[id].annotations
A.coupling_annotations(m::Model, id::String) = m.couplings[id].annotations
A.reaction_notes(m::Model, id::String) = m.reactions[id].notes
A.metabolite_notes(m::Model, id::String) = m.metabolites[id].notes
A.gene_notes(m::Model, id::String) = m.genes[id].notes
A.coupling_notes(m::Model, id::String) = m.couplings[id].notes

gene_symbol(m::Model, id::String) = m.genes[id].symbol

function A.stoichiometry(m::Model)
    midxs = Dict(mid => idx for (idx, mid) in enumerate(A.metabolites(m)))
    I = Int[]
    J = Int[]
    V = Float64[]
    for (ridx, rid) in enumerate(A.reactions(m))
        for (smid, v) in m.reactions[rid].stoichiometry
            push!(I, midxs[smid])
            push!(J, ridx)
            push!(V, v)
        end
    end
    S.sparse(I, J, V, A.n_metabolites(m), A.n_reactions(m))
end

function A.coupling(m::Model)
    ridxs = Dict(rid => idx for (idx, rid) in enumerate(A.reactions(m)))
    I = Int[]
    J = Int[]
    V = Float64[]
    for (cidx, cid) in enumerate(A.couplings(m))
        for (rid, v) in m.couplings[cid].reaction_weights
            push!(I, cidx)
            push!(J, ridxs[rid])
            push!(V, v)
        end
    end
    S.sparse(I, J, V, A.n_couplings(m), A.n_reactions(m))
end

A.bounds(m::Model) = (
    [m.reactions[rid].lower_bound for rid in A.reactions(m)],
    [m.reactions[rid].upper_bound for rid in A.reactions(m)],
)

A.coupling_bounds(m::Model) = (
    [m.couplings[cid].lower_bound for cid in A.couplings(m)],
    [m.couplings[cid].upper_bound for cid in A.couplings(m)],
)

A.balance(m::Model) =
    S.sparse(Float64[m.metabolites[mid].balance for mid in A.metabolites(m)])
A.objective(m::Model) =
    S.sparse(Float64[m.reactions[rid].objective_coefficient for rid in A.reactions(m)])

A.reaction_gene_association_dnf(m::Model, id::String) = begin
    grrs = m.reactions[id].gene_association
    isnothing(grrs) && return nothing
    [collect(keys(iso.gene_product_stoichiometry)) for iso in grrs]
end

A.reaction_gene_products_available(m::Model, id::String, fn::Function) =
    A.reaction_gene_products_available_from_dnf(m, id, fn) # TODO

A.reaction_stoichiometry(m::Model, id::String) = m.reactions[id].stoichiometry

A.metabolite_formula(m::Model, id::String) = m.metabolites[id].formula
A.metabolite_charge(m::Model, id::String) = m.metabolites[id].charge
A.metabolite_compartment(m::Model, id::String) = m.metabolites[id].compartment

A.coupling_weights(m::Model, id::String) = m.couplings[id].reaction_weights

A.load(::Type{Model}, path::String) = S.deserialize(path)
A.save(m::Model, path::String) = S.serialize(path, m)
A.filename_extensions(::Type{Model}) = ["vibrio-natriegens-serialized-fbc"]

function Base.convert(::Type{Model}, x::A.AbstractFBCModel)
    (lbs, ubs) = A.bounds(x)
    (clbs, cubs) = A.coupling_bounds(x)
    Model(
        reactions = Dict(
            r => Reaction(
                name = A.reaction_name(x, r),
                lower_bound = lb,
                upper_bound = ub,
                stoichiometry = A.reaction_stoichiometry(x, r),
                objective_coefficient = o,
                annotations = A.reaction_annotations(x, r),
                notes = A.reaction_notes(x, r),
            ) for (r, o, lb, ub) in zip(A.reactions(x), A.objective(x), lbs, ubs)
        ),
        metabolites = Dict(
            m => Metabolite(
                name = A.metabolite_name(x, m),
                balance = b,
                formula = A.metabolite_formula(x, m),
                charge = A.metabolite_charge(x, m),
                compartment = A.metabolite_compartment(x, m),
                annotations = A.metabolite_annotations(x, m),
                notes = A.metabolite_notes(x, m),
            ) for (m, b) in zip(A.metabolites(x), A.balance(x))
        ),
        genes = Dict(
            g => Gene(
                name = A.gene_name(x, g),
                annotations = A.gene_annotations(x, g),
                notes = A.gene_notes(x, g),
            ) for g in A.genes(x)
        ),
        couplings = Dict(
            c => Coupling(
                name = A.coupling_name(x, c),
                lower_bound = lb,
                upper_bound = ub,
                reaction_weights = A.coupling_weights(x, c),
                annotations = A.coupling_annotations(x, c),
                notes = A.coupling_notes(x, c),
            ) for (c, lb, ub) in zip(A.couplings(x), clbs, cubs)
        ),
    )
end
